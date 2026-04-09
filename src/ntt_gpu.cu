#include "ntt_gpu.cuh"
#include <cuda_runtime.h>
#include <iostream>
#include <cassert>

#define CHECK_CUDA(call) assert((call) == cudaSuccess)

// -----------------------------------------------
//              Helper Functions
// -----------------------------------------------
__device__ inline uint32_t mod_add_gpu(uint32_t a, uint32_t b, uint32_t q) {
    uint32_t sum = a + b;
    return (sum >= q) ? (sum - q) : sum;
}

__device__ inline uint32_t mod_sub_gpu(uint32_t a, uint32_t b, uint32_t q) {
    return (a >= b) ? (a - b) : (a + q - b);
}

template <bool UseBarrett>
__device__ inline uint32_t mod_mul_gpu(uint32_t a, uint32_t b, uint32_t q, uint64_t mu) {
    uint64_t x = static_cast<uint64_t>(a) * b;
    if constexpr (UseBarrett) {
        uint64_t q_est = __umul64hi(x, mu);
        uint32_t r = static_cast<uint32_t>(x - q_est * q);
        return (r >= q) ? (r - q) : r;
    } else {
        return static_cast<uint32_t>(x % q);
    }
}

// -----------------------------------------------
//              Radix 2
// -----------------------------------------------

// stage independent Gentleman Sande
template <bool UseBarrett>
__global__ void dif_gs_ntt_stage(uint32_t* a, const uint32_t* twiddles, uint32_t N, uint32_t q, uint64_t mu, uint32_t len) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N / 2) return;

    // compute indicies
    uint32_t half_len = len / 2;
    uint32_t group = tid / half_len;
    uint32_t j = tid % half_len;
    uint32_t i = group * len;

    // compute stride
    uint32_t step = N / len;

    // fetch precomputed twiddleu
    uint32_t w = twiddles[j * step];

    uint32_t u = a[i + j];
    uint32_t v = a[i + j + half_len];

    // butterfly
    a[i + j] = mod_add_gpu(u, v, q);
    uint32_t diff = mod_sub_gpu(u, v, q);
    a[i + j + half_len] = mod_mul_gpu<UseBarrett>(diff, w, q, mu);
}

// stage independent Cooley Tukey
template <bool UseBarrett>
__global__ void dit_ct_intt_stage(uint32_t* a, const uint32_t* inv_twiddles, uint32_t N, uint32_t q, uint64_t mu, uint32_t len) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N / 2) return;

    uint32_t half_len = len / 2;
    uint32_t group = tid / half_len;
    uint32_t j = tid % half_len;
    uint32_t i = group * len;
    uint32_t step = N / len;
    uint32_t w = inv_twiddles[j * step];

    uint32_t u = a[i + j];
    uint32_t v = mod_mul_gpu<UseBarrett>(a[i + j + half_len], w, q, mu);

    a[i + j] = mod_add_gpu(u, v, q);
    a[i + j + half_len] = mod_sub_gpu(u, v, q);
}

template <bool UseBarrett>
__global__ void normalize_intt(uint32_t* a, uint32_t N, uint32_t q, uint64_t mu, uint32_t inv_N) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < N) {
        a[tid] = mod_mul_gpu<UseBarrett>(a[tid], inv_N, q, mu);
    }
}

// -----------------------------------------------
//              Radix 4
// -----------------------------------------------

template <bool UseBarrett>
__global__ void dif_radix4_ntt_stage(uint32_t* a, const uint32_t* twiddles, uint32_t N, uint32_t q, uint64_t mu, uint32_t len, uint32_t W_4) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N / 4) return;

    uint32_t quarter_len = len / 4;
    uint32_t group = tid / quarter_len;
    uint32_t j = tid % quarter_len;
    uint32_t i = group * len;
    uint32_t step = N / len;
    uint32_t W1_idx = j * step;
    
    uint32_t w1 = twiddles[W1_idx];
    uint32_t w2 = twiddles[W1_idx * 2];
    uint32_t w3 = twiddles[W1_idx * 3];

    uint32_t u0 = a[i + j];
    uint32_t u1 = a[i + j + quarter_len];
    uint32_t u2 = a[i + j + 2 * quarter_len];
    uint32_t u3 = a[i + j + 3 * quarter_len];

    uint32_t A = mod_add_gpu(u0, u2, q);
    uint32_t B = mod_sub_gpu(u0, u2, q);
    uint32_t C = mod_add_gpu(u1, u3, q);
    uint32_t D = mod_sub_gpu(u1, u3, q);

    D = mod_mul_gpu<UseBarrett>(D, W_4, q, mu);

    uint32_t out0 = mod_add_gpu(A, C, q);
    uint32_t out1 = mod_add_gpu(B, D, q);
    uint32_t out2 = mod_sub_gpu(A, C, q);
    uint32_t out3 = mod_sub_gpu(B, D, q);

    a[i + j]                   = out0;
    a[i + j + quarter_len]     = mod_mul_gpu<UseBarrett>(out1, w1, q, mu);
    a[i + j + 2 * quarter_len] = mod_mul_gpu<UseBarrett>(out2, w2, q, mu);
    a[i + j + 3 * quarter_len] = mod_mul_gpu<UseBarrett>(out3, w3, q, mu);
}

template <bool UseBarrett>
__global__ void dit_radix4_intt_stage(uint32_t* a, const uint32_t* inv_twiddles, uint32_t N, uint32_t q, uint64_t mu, uint32_t len, uint32_t inv_W_4) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N / 4) return;

    uint32_t quarter_len = len / 4;
    uint32_t group = tid / quarter_len;
    uint32_t j = tid % quarter_len;
    uint32_t i = group * len;
    uint32_t step = N / len;
    uint32_t W1_idx = j * step;
    
    uint32_t w1 = inv_twiddles[W1_idx];
    uint32_t w2 = inv_twiddles[W1_idx * 2];
    uint32_t w3 = inv_twiddles[W1_idx * 3];

    uint32_t u0 = a[i + j];
    uint32_t u1 = mod_mul_gpu<UseBarrett>(a[i + j + quarter_len], w1, q, mu);
    uint32_t u2 = mod_mul_gpu<UseBarrett>(a[i + j + 2 * quarter_len], w2, q, mu);
    uint32_t u3 = mod_mul_gpu<UseBarrett>(a[i + j + 3 * quarter_len], w3, q, mu);

    uint32_t A = mod_add_gpu(u0, u2, q);
    uint32_t B = mod_sub_gpu(u0, u2, q);
    uint32_t C = mod_add_gpu(u1, u3, q);
    uint32_t D = mod_sub_gpu(u1, u3, q);

    D = mod_mul_gpu<UseBarrett>(D, inv_W_4, q, mu);

    a[i + j]                   = mod_add_gpu(A, C, q);
    a[i + j + quarter_len]     = mod_add_gpu(B, D, q);
    a[i + j + 2 * quarter_len] = mod_sub_gpu(A, C, q);
    a[i + j + 3 * quarter_len] = mod_sub_gpu(B, D, q);
}



// ----------------------------
//    Host side kernel launch
// ----------------------------

float ntt_gpu_dif(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t root_of_unity, bool use_barrett) {
    uint32_t N = data.size();
    uint32_t bytes = N * sizeof(uint32_t);
    uint32_t *d_data, *d_twiddles;

    CHECK_CUDA(cudaMalloc(&d_data, bytes));
    CHECK_CUDA(cudaMalloc(&d_twiddles, bytes));
    
    // precompute twiddles
    std::vector<uint32_t> twiddles(N);
    uint64_t current_W = 1;
    for (uint32_t i = 0; i < N; i++) {
        twiddles[i] = static_cast<uint32_t>(current_W);
        current_W = (current_W * root_of_unity) % q;
    }

    CHECK_CUDA(cudaMemcpy(d_data, data.data(), bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_twiddles, twiddles.data(), bytes, cudaMemcpyHostToDevice));

    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));

    int threads = 256;
    int blocks = (N / 2 + threads - 1) / threads;

    CHECK_CUDA(cudaEventRecord(start));
    for (uint32_t len = N; len >= 2; len >>= 1) {
        if (use_barrett) dif_gs_ntt_stage<true><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
        else             dif_gs_ntt_stage<false><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
    }
    CHECK_CUDA(cudaEventRecord(stop));
    CHECK_CUDA(cudaEventSynchronize(stop));

    float ms = 0;
    CHECK_CUDA(cudaEventElapsedTime(&ms, start, stop));

    CHECK_CUDA(cudaMemcpy(data.data(), d_data, bytes, cudaMemcpyDeviceToHost));
    cudaFree(d_data); cudaFree(d_twiddles);
    cudaEventDestroy(start); cudaEventDestroy(stop);
    return ms;
}

float intt_gpu_dit(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t inv_root_of_unity, uint32_t inv_N, bool use_barrett) {
    uint32_t N = data.size();
    uint32_t bytes = N * sizeof(uint32_t);
    uint32_t *d_data, *d_twiddles;

    CHECK_CUDA(cudaMalloc(&d_data, bytes));
    CHECK_CUDA(cudaMalloc(&d_twiddles, bytes));

    std::vector<uint32_t> twiddles(N);
    uint64_t current_W = 1;
    for (uint32_t i = 0; i < N; i++) {
        twiddles[i] = static_cast<uint32_t>(current_W);
        current_W = (current_W * inv_root_of_unity) % q;
    }

    CHECK_CUDA(cudaMemcpy(d_data, data.data(), bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_twiddles, twiddles.data(), bytes, cudaMemcpyHostToDevice));

    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));

    int threads = 256;
    int blocks = (N / 2 + threads - 1) / threads;

    CHECK_CUDA(cudaEventRecord(start));
    for (uint32_t len = 2; len <= N; len <<= 1) {
        if (use_barrett) dit_ct_intt_stage<true><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
        else             dit_ct_intt_stage<false><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
    }
    
    int normBlocks = (N + threads - 1) / threads;
    if (use_barrett) normalize_intt<true><<<normBlocks, threads>>>(d_data, N, q, mu, inv_N);
    else             normalize_intt<false><<<normBlocks, threads>>>(d_data, N, q, mu, inv_N);
    
    CHECK_CUDA(cudaEventRecord(stop));
    CHECK_CUDA(cudaEventSynchronize(stop));

    float ms = 0;
    CHECK_CUDA(cudaEventElapsedTime(&ms, start, stop));

    CHECK_CUDA(cudaMemcpy(data.data(), d_data, bytes, cudaMemcpyDeviceToHost));
    cudaFree(d_data); cudaFree(d_twiddles);
    cudaEventDestroy(start); cudaEventDestroy(stop);
    return ms;
}

float ntt_gpu_radix4_dif(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t root_of_unity, bool use_barrett) {
    uint32_t N = data.size();
    uint32_t bytes = N * sizeof(uint32_t);
    uint32_t *d_data, *d_twiddles;

    CHECK_CUDA(cudaMalloc(&d_data, bytes));
    CHECK_CUDA(cudaMalloc(&d_twiddles, bytes));
    
    std::vector<uint32_t> twiddles(N);
    uint64_t current_W = 1;
    for (uint32_t i = 0; i < N; i++) {
        twiddles[i] = static_cast<uint32_t>(current_W);
        current_W = (current_W * root_of_unity) % q;
    }
    CHECK_CUDA(cudaMemcpy(d_data, data.data(), bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_twiddles, twiddles.data(), bytes, cudaMemcpyHostToDevice));

    uint32_t W_4 = twiddles[N / 4];
    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));

    int threads = 256;
    uint32_t len = N;

    CHECK_CUDA(cudaEventRecord(start));
    // run radix 4 
    while (len >= 4) {
        int blocks = (N / 4 + threads - 1) / threads;
        if (use_barrett) dif_radix4_ntt_stage<true><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len, W_4);
        else             dif_radix4_ntt_stage<false><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len, W_4);
        len >>= 2;
    }
    // fallback to Radix 2 if N is an even power of 2
    if (len == 2) {
        int blocks = (N / 2 + threads - 1) / threads;
        if (use_barrett) dif_gs_ntt_stage<true><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
        else             dif_gs_ntt_stage<false><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
    }
    CHECK_CUDA(cudaEventRecord(stop));
    CHECK_CUDA(cudaEventSynchronize(stop));

    float ms = 0;
    CHECK_CUDA(cudaEventElapsedTime(&ms, start, stop));

    CHECK_CUDA(cudaMemcpy(data.data(), d_data, bytes, cudaMemcpyDeviceToHost));
    cudaFree(d_data); cudaFree(d_twiddles);
    cudaEventDestroy(start); cudaEventDestroy(stop);
    return ms;
}

float intt_gpu_radix4_dit(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t inv_root_of_unity, uint32_t inv_N, bool use_barrett) {
    uint32_t N = data.size();
    uint32_t bytes = N * sizeof(uint32_t);
    uint32_t *d_data, *d_twiddles;

    CHECK_CUDA(cudaMalloc(&d_data, bytes));
    CHECK_CUDA(cudaMalloc(&d_twiddles, bytes));
    
    std::vector<uint32_t> twiddles(N);
    uint64_t current_W = 1;
    for (uint32_t i = 0; i < N; i++) {
        twiddles[i] = static_cast<uint32_t>(current_W);
        current_W = (current_W * inv_root_of_unity) % q;
    }
    CHECK_CUDA(cudaMemcpy(d_data, data.data(), bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_twiddles, twiddles.data(), bytes, cudaMemcpyHostToDevice));

    uint32_t inv_W_4 = twiddles[N / 4]; 
    int stages = 0; uint32_t temp = N;
    while (temp > 1) { temp >>= 1; stages++; }

    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));

    int threads = 256;
    uint32_t len;

    CHECK_CUDA(cudaEventRecord(start));
    // if N is an odd power of 2, stages is odd
    // process the first stage using radix 2
    if (stages % 2 != 0) {
        len = 2;
        int blocks = (N / 2 + threads - 1) / threads;
        if (use_barrett) dit_ct_intt_stage<true><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);
        else             dit_ct_intt_stage<false><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len);

        // next stage length 4 * 2 = 8
        len = 8; 
    } else {
        len = 4;
    }

    // continue processing using radix 4
    while (len <= N) {
        int blocks = (N / 4 + threads - 1) / threads;
        if (use_barrett) dit_radix4_intt_stage<true><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len, inv_W_4);
        else             dit_radix4_intt_stage<false><<<blocks, threads>>>(d_data, d_twiddles, N, q, mu, len, inv_W_4);
        len <<= 2; 
    }

    int normBlocks = (N + threads - 1) / threads;
    if (use_barrett) normalize_intt<true><<<normBlocks, threads>>>(d_data, N, q, mu, inv_N);
    else             normalize_intt<false><<<normBlocks, threads>>>(d_data, N, q, mu, inv_N);
    
    CHECK_CUDA(cudaEventRecord(stop));
    CHECK_CUDA(cudaEventSynchronize(stop));

    float ms = 0;
    CHECK_CUDA(cudaEventElapsedTime(&ms, start, stop));

    CHECK_CUDA(cudaMemcpy(data.data(), d_data, bytes, cudaMemcpyDeviceToHost));
    cudaFree(d_data); cudaFree(d_twiddles);
    cudaEventDestroy(start); cudaEventDestroy(stop);
    return ms;
}