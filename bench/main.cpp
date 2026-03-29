#include <openfhe.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>

#include "utils.h"
#include "ntt_cpu.h"
#include "ntt_openfhe.h"

#define NUM_RUNS 5

using namespace lbcrypto;
using namespace std;

int main() {
  uint32_t N = 1 << 10;
  uint32_t cyclotomic_order = 2 * N;
  uint32_t bit_size = 28;
  uint32_t num_limbs = 128;

  printf("| N = %u | L = %u |\n\n", N, num_limbs);

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

  double naive_time = 0;
  double fast_time = 0;
  double prod_time = 0;

  for (uint32_t run = 0; run < NUM_RUNS; run++) {
    for (uint32_t i = 0; i < num_limbs; i++) {
      std::vector<uint32_t> poly = original_rns_poly[i];
      uint32_t q = rns_params[i].q;

      // generate twiddles
      std::vector<uint32_t> W = generate_sequential_twiddles(N, q, rns_params[i].root);
      std::vector<uint32_t> inv_W = generate_sequential_twiddles(N, q, rns_params[i].inv_root);

      // naive
      auto start = std::chrono::high_resolution_clock::now();
      std::vector<uint32_t> res_naive = naive_ntt(poly, q, rns_params[i].mu, W);
      res_naive = naive_intt(res_naive, q, rns_params[i].mu, inv_W, rns_params[i].inv_N);
      auto end = std::chrono::high_resolution_clock::now();
      naive_time += std::chrono::duration<double, std::milli>(end - start).count();

      // fast
      start = std::chrono::high_resolution_clock::now();
      std::vector<uint32_t> res_fast = fast_gs_ntt(poly, q, rns_params[i].mu, rns_params[i].root);
      res_fast = fast_ct_intt(res_fast, q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N);
      end = std::chrono::high_resolution_clock::now();
      fast_time += std::chrono::duration<double, std::milli>(end - start).count();

      // production
      auto start_prod = std::chrono::high_resolution_clock::now();
      std::vector<uint32_t> res_prod = prod_gs_ntt(poly, q, rns_params[i].mu, rns_params[i].omega_pow);
      res_prod = prod_ct_intt(res_prod, q, rns_params[i].mu, rns_params[i].inv_omega_pow, rns_params[i].inv_N);
      auto end_prod = std::chrono::high_resolution_clock::now();
      prod_time += std::chrono::duration<double, std::milli>(end_prod - start_prod).count();

      // correctness check for custom CPU
      for (uint32_t j = 0; j < N; j++) {
        if (res_naive[j] != original_rns_poly[i][j]) {
          printf("Naive NTT Error, mismatch at limb %d idx %d \n", i, j);
          return 1;
        }
        if (res_fast[j] != original_rns_poly[i][j]) {
          printf("Fast NTT Error, mismatch at limb %d idx %d \n", i, j);
          return 1;
        }
        if (res_prod[j] != original_rns_poly[i][j]) {
          printf("Prod NTT Error, mismatch at limb %d idx %d \n", i, j);
          return 1;
        }
      }
    }
  }

  // run OpenFHE + OMP and check correctness
  OpenFHEBenchResult openfhe_res = benchmark_openfhe_rns_ntt(original_rns_poly, N, cyclotomic_order, rns_params, NUM_RUNS);
  double openfhe_time = openfhe_res.time_ms;

  for (uint32_t i = 0; i < num_limbs; i++) {
    for (uint32_t j = 0; j < N; j++) {
      if (openfhe_res.ntt_result[i][j] != original_rns_poly[i][j]) {
        printf("OpenFHE NTT Error, mismatch at limb %d idx %d \n", i, j);
        return 1;
      }
    }
  }

  naive_time /= NUM_RUNS;
  fast_time /= NUM_RUNS;
  prod_time /= NUM_RUNS;

  printf("Average Naive NTT Time over %d runs: %f ms\n", NUM_RUNS, naive_time);
  printf("Average Fast NTT Time over %d runs: %f ms\n", NUM_RUNS, fast_time);
  printf("Average Prod NTT Time over %d runs: %f ms\n", NUM_RUNS, prod_time);
  printf("Average OpenFHE NTT Time over %d runs: %f ms\n", NUM_RUNS, openfhe_time);

  return 0;
}