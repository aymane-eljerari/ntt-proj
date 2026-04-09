#include "ntt_cpu.h"
#include "utils.h"

using namespace std;

/* ----------------------
    naive NTT O(N^2) 
------------------------- */
template <bool UseBarrett>
vector<uint32_t> naive_ntt(vector<uint32_t> a, uint32_t q, uint64_t mu, const vector<uint32_t>& W) {
    uint32_t N = a.size();
    vector<uint32_t> result(N, 0);

    for (uint32_t i = 0; i < N; i++) {
        uint32_t sum = 0;
        for (uint32_t j = 0; j < N; j++) {
            uint32_t idx = (i * j) % N;
            uint32_t coef = mod_mul<UseBarrett>(a[j], W[idx], q, mu);
            sum = mod_add<UseBarrett>(sum, coef, q);
        }
        result[i] = sum;
    }
    return result;
}

template <bool UseBarrett>
vector<uint32_t> naive_intt(vector<uint32_t> a, uint32_t q, uint64_t mu, const vector<uint32_t>& inv_W, const uint32_t inv_N) {
    a = naive_ntt<UseBarrett>(a, q, mu, inv_W);
    for (uint32_t i = 0; i < a.size(); i++) {
        a[i] = mod_mul<UseBarrett>(a[i], inv_N, q, mu);
    }
    return a;
}

/*
------------------------
    fast NTT O(n log n)
------------------------
*/
template <bool UseBarrett>
vector<uint32_t> fast_gs_ntt(vector<uint32_t> a, uint32_t q, uint64_t mu, uint32_t root) {
    uint32_t N = a.size();
    for (uint32_t len = N; len >= 2; len >>= 1) {
        uint32_t wlen = mod_exp<UseBarrett>(root, N / len, q, mu);
        for (uint32_t i = 0; i < N; i += len) {
            uint32_t w = 1;
            for (uint32_t j = 0; j < len / 2; j++) {
                uint32_t u = a[i + j];
                uint32_t v = a[i + j + len / 2];
                a[i + j] = mod_add<UseBarrett>(u, v, q);
                uint32_t diff = mod_sub<UseBarrett>(u, v, q);
                
                a[i + j + len / 2] = mod_mul<UseBarrett>(diff, w, q, mu);
                w = mod_mul<UseBarrett>(w, wlen, q, mu);
            }
        }
    }
    return a;
}

template <bool UseBarrett>
vector<uint32_t> fast_ct_intt(vector<uint32_t> a, uint32_t q, uint64_t mu, uint32_t inv_root, uint32_t inv_N) {
    uint32_t N = a.size();
    for (uint32_t len = 2; len <= N; len <<= 1) {
        uint32_t wlen = mod_exp<UseBarrett>(inv_root, N / len, q, mu);
        for (uint32_t i = 0; i < N; i += len) {
            uint32_t w = 1;
            for (uint32_t j = 0; j < len / 2; j++) {
                uint32_t u = a[i + j];
                uint32_t v = mod_mul<UseBarrett>(a[i + j + len / 2], w, q, mu);
                
                a[i + j] = mod_add<UseBarrett>(u, v, q);
                a[i + j + len / 2] = mod_sub<UseBarrett>(u, v, q);

                w = mod_mul<UseBarrett>(w, wlen, q, mu);
            }
        }
    }
    for (uint32_t i = 0; i < N; i++) {
        a[i] = mod_mul<UseBarrett>(a[i], inv_N, q, mu);
    }
    return a;
}

/*
-----------------------------------------
    production NTT (precomputed twiddles)
-----------------------------------------
*/
template <bool UseBarrett>
vector<uint32_t> prod_gs_ntt(vector<uint32_t> a, uint32_t q, uint64_t mu, const vector<uint32_t>& omega_pow) {
    uint32_t N = a.size();
    for (uint32_t len = N; len >= 2; len >>= 1) {
        uint32_t step = N / len;
        uint32_t half_len = len / 2;
        for (uint32_t i = 0; i < N; i += len) {
            for (uint32_t j = 0; j < half_len; j++) {
                uint32_t w = omega_pow[j * step];
                
                uint32_t u = a[i + j];
                uint32_t v = a[i + j + half_len];
                
                a[i + j] = mod_add<UseBarrett>(u, v, q);
                uint32_t diff = mod_sub<UseBarrett>(u, v, q);
                a[i + j + half_len] = mod_mul<UseBarrett>(diff, w, q, mu);
            }
        }
    }
    return a;
}

template <bool UseBarrett>
vector<uint32_t> prod_ct_intt(vector<uint32_t> a, uint32_t q, uint64_t mu, const vector<uint32_t>& inv_omega_pow, uint32_t inv_N) {
    uint32_t N = a.size();
    for (uint32_t len = 2; len <= N; len <<= 1) {
        uint32_t step = N / len;
        uint32_t half_len = len / 2;
        for (uint32_t i = 0; i < N; i += len) {
            for (uint32_t j = 0; j < half_len; j++) {
                uint32_t w = inv_omega_pow[j * step];
                
                uint32_t u = a[i + j];
                uint32_t v = mod_mul<UseBarrett>(a[i + j + half_len], w, q, mu);
                
                a[i + j] = mod_add<UseBarrett>(u, v, q);
                a[i + j + half_len] = mod_sub<UseBarrett>(u, v, q);
            }
        }
    }
    
    for (uint32_t i = 0; i < N; i++) {
        a[i] = mod_mul<UseBarrett>(a[i], inv_N, q, mu);
    }
    return a;
}

// instanciate templates (fixes linker issue)
template vector<uint32_t> naive_ntt<true>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&);
template vector<uint32_t> naive_ntt<false>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&);

template vector<uint32_t> naive_intt<true>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&, uint32_t);
template vector<uint32_t> naive_intt<false>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&, uint32_t);

template vector<uint32_t> fast_gs_ntt<true>(vector<uint32_t>, uint32_t, uint64_t, uint32_t);
template vector<uint32_t> fast_gs_ntt<false>(vector<uint32_t>, uint32_t, uint64_t, uint32_t);

template vector<uint32_t> fast_ct_intt<true>(vector<uint32_t>, uint32_t, uint64_t, uint32_t, uint32_t);
template vector<uint32_t> fast_ct_intt<false>(vector<uint32_t>, uint32_t, uint64_t, uint32_t, uint32_t);

template vector<uint32_t> prod_gs_ntt<true>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&);
template vector<uint32_t> prod_gs_ntt<false>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&);

template vector<uint32_t> prod_ct_intt<true>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&, uint32_t);
template vector<uint32_t> prod_ct_intt<false>(vector<uint32_t>, uint32_t, uint64_t, const vector<uint32_t>&, uint32_t);