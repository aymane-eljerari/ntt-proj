#include <openfhe.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <omp.h> 

#include "utils.h"
#include "ntt_cpu.h"
#include "ntt_openfhe.h"
#include "ntt_gpu.cuh"
#include "arg_parser.h"

#define NUM_RUNS 2

using namespace lbcrypto;
using namespace std;

// Struct to hold and return benchmark results
struct BenchmarkResults {
  double naive_time_ms;
  double fast_time_ms;
  double prod_time_ms;
  double gpu_radix2_time_ms;
  double gpu_radix4_time_ms;
  double openfhe_time_ms;
};

int main(int argc, char** argv) {
  BenchConfig config = parse_args(argc, argv);

  uint32_t N = config.N;
  uint32_t cyclotomic_order = 2 * N;
  uint32_t bit_size = 28;
  uint32_t num_limbs = config.num_limbs;

  // Set OpenMP threads
  omp_set_num_threads(config.threads);

  printf("| N = %u | L = %u | Threads = %u |\n\n", N, num_limbs, config.threads);

  // init rns polynomial + get the first prime
  std::vector<RNSLimbParams> rns_params(num_limbs);
  NativeInteger current_prime = FirstPrime<NativeInteger>(bit_size, cyclotomic_order);

  // generate and store rns parameters
  for (uint32_t i = 0; i < num_limbs; i++) {
    NativeInteger psi_unity = RootOfUnity<NativeInteger>(cyclotomic_order, current_prime);
    uint32_t q = (uint32_t)current_prime.ConvertToInt();
    uint32_t omega = mod_exp((uint32_t)psi_unity.ConvertToInt(), 2, q);
    uint32_t inv_omega = mod_inverse(omega, q);
    uint32_t inv_N = mod_inverse(N, q);

    rns_params[i].q = q;
    rns_params[i].mu = (uint64_t)(((__uint128_t) 1 << 64) / q);
    rns_params[i].root = omega;
    rns_params[i].inv_root = inv_omega;
    rns_params[i].inv_N = inv_N;

    rns_params[i].omega_pow = generate_sequential_twiddles(N, q, omega);
    rns_params[i].inv_omega_pow = generate_sequential_twiddles(N, q, inv_omega);
    
    current_prime = NextPrime<NativeInteger>(current_prime, cyclotomic_order);
  }

  // rng sampling to populate the rns limbs
  mt19937 gen(42);
  std::vector<std::vector<uint32_t>> original_rns_poly(num_limbs, vector<uint32_t>(N));
  for (uint32_t i = 0; i < num_limbs; i++) {
    uniform_int_distribution<uint32_t> dist(0, rns_params[i].q - 1);
    for (uint32_t j = 0; j < N; j++) {
        original_rns_poly[i][j] = dist(gen);
    }
  }

  BenchmarkResults results = {0};
  bool run_all = config.run_all();

  for (uint32_t run = 0; run < NUM_RUNS; run++) {
    for (uint32_t i = 0; i < num_limbs; i++) {
      std::vector<uint32_t> poly = original_rns_poly[i];
      uint32_t q = rns_params[i].q;

      // generate Twiddles
      std::vector<uint32_t> W = generate_sequential_twiddles(N, q, rns_params[i].root);
      std::vector<uint32_t> inv_W = generate_sequential_twiddles(N, q, rns_params[i].inv_root);

      // Naive
      if (run_all || config.run_cpu_naive) {
          auto start = std::chrono::high_resolution_clock::now();
          std::vector<uint32_t> res_naive = naive_ntt(poly, q, rns_params[i].mu, W);
          res_naive = naive_intt(res_naive, q, rns_params[i].mu, inv_W, rns_params[i].inv_N);
          auto end = std::chrono::high_resolution_clock::now();
          results.naive_time_ms += std::chrono::duration<double, std::milli>(end - start).count();

          // Correctness Check
          for (uint32_t j = 0; j < N; j++) {
            if (res_naive[j] != original_rns_poly[i][j]) {
              printf("Naive NTT Error, mismatch at limb %d idx %d \n", i, j);
              return 1;
            }
          }
      }

      // fast
      if (run_all || config.run_cpu_fast) {
          auto start = std::chrono::high_resolution_clock::now();
          std::vector<uint32_t> res_fast = fast_gs_ntt(poly, q, rns_params[i].mu, rns_params[i].root);
          res_fast = fast_ct_intt(res_fast, q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N);
          auto end = std::chrono::high_resolution_clock::now();
          results.fast_time_ms += std::chrono::duration<double, std::milli>(end - start).count();

          for (uint32_t j = 0; j < N; j++) {
            if (res_fast[j] != original_rns_poly[i][j]) {
              printf("Fast NTT Error, mismatch at limb %d idx %d \n", i, j);
              return 1;
            }
          }
      }

      // Production
      if (run_all || config.run_cpu_prod) {
          auto start = std::chrono::high_resolution_clock::now();
          std::vector<uint32_t> res_prod = prod_gs_ntt(poly, q, rns_params[i].mu, rns_params[i].omega_pow);
          res_prod = prod_ct_intt(res_prod, q, rns_params[i].mu, rns_params[i].inv_omega_pow, rns_params[i].inv_N);
          auto end = std::chrono::high_resolution_clock::now();
          results.prod_time_ms += std::chrono::duration<double, std::milli>(end - start).count();

          for (uint32_t j = 0; j < N; j++) {
            if (res_prod[j] != original_rns_poly[i][j]) {
              printf("Prod NTT Error, mismatch at limb %d idx %d \n", i, j);
              return 1;
            }
          }
      }

      // GPU radix 2
      if (run_all || config.run_gpu_radix2) {
          std::vector<uint32_t> res_gpu_radix2 = poly;
          auto start = std::chrono::high_resolution_clock::now();
          ntt_gpu_dif(res_gpu_radix2, q, rns_params[i].mu, rns_params[i].root);      
          intt_gpu_dit(res_gpu_radix2, q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N); 
          auto end = std::chrono::high_resolution_clock::now();
          results.gpu_radix2_time_ms += std::chrono::duration<double, std::milli>(end - start).count();

          for (uint32_t j = 0; j < N; j++) {
            if (res_gpu_radix2[j] != original_rns_poly[i][j]) {
              printf("GPU Radix-2 NTT Error, mismatch at limb %d idx %d \n", i, j);
              return 1;
            }
          }
      }

      // GPU radix 4
      if (run_all || config.run_gpu_radix4) {
          std::vector<uint32_t> res_gpu_radix4 = poly;
          auto start = std::chrono::high_resolution_clock::now();
          ntt_gpu_radix4_dif(res_gpu_radix4, q, rns_params[i].mu, rns_params[i].root);      
          intt_gpu_radix4_dit(res_gpu_radix4, q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N); 
          auto end = std::chrono::high_resolution_clock::now();
          results.gpu_radix4_time_ms += std::chrono::duration<double, std::milli>(end - start).count();

          for (uint32_t j = 0; j < N; j++) {
            if (res_gpu_radix4[j] != original_rns_poly[i][j]) {
              printf("GPU Radix-4 NTT Error, mismatch at limb %d idx %d \n", i, j);
              return 1;
            }
          }
      }
    }
  }

  // OpenFHE always runs to ensure baseline comparison
  OpenFHEBenchResult openfhe_res = benchmark_openfhe_rns_ntt(original_rns_poly, N, cyclotomic_order, rns_params, NUM_RUNS);
  results.openfhe_time_ms = openfhe_res.time_ms;

  for (uint32_t i = 0; i < num_limbs; i++) {
    for (uint32_t j = 0; j < N; j++) {
      if (openfhe_res.ntt_result[i][j] != original_rns_poly[i][j]) {
        printf("OpenFHE NTT Error, mismatch at limb %d idx %d \n", i, j);
        return 1;
      }
    }
  }

  results.naive_time_ms /= NUM_RUNS;
  results.fast_time_ms /= NUM_RUNS;
  results.prod_time_ms /= NUM_RUNS;
  results.gpu_radix2_time_ms /= NUM_RUNS;
  results.gpu_radix4_time_ms /= NUM_RUNS;

  printf("================ BENCHMARK RESULTS ================\n");
  if (run_all || config.run_cpu_naive) printf("Average Naive NTT Time      : %f ms\n", results.naive_time_ms);
  if (run_all || config.run_cpu_fast)  printf("Average Fast NTT Time       : %f ms\n", results.fast_time_ms);
  if (run_all || config.run_cpu_prod)  printf("Average Prod NTT Time       : %f ms\n", results.prod_time_ms);
  if (run_all || config.run_gpu_radix2)printf("Average GPU Radix-2 Time    : %f ms\n", results.gpu_radix2_time_ms);
  if (run_all || config.run_gpu_radix4)printf("Average GPU Radix-4 Time    : %f ms\n", results.gpu_radix4_time_ms);
  printf("Average OpenFHE NTT Time    : %f ms\n", results.openfhe_time_ms);
  printf("===================================================\n");

  return 0;
}