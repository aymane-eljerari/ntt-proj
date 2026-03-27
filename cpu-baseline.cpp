#include <openfhe.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>

using namespace lbcrypto;
using namespace std;

struct RNSLimbParams {
  uint32_t q;
  uint32_t inv_N;
  uint32_t root;
  uint32_t inv_root;
  vector<uint32_t> omega_pow;
  vector<uint32_t> inv_omega_pow;
};


/* ----------------------
    helper functions 
------------------------- */

uint32_t mod_exp(uint32_t base, uint32_t exp, uint32_t mod) {
  uint32_t res = 1;
  base = base % mod;
  while (exp > 0) {
    if (exp % 2 == 1) res = (uint32_t)(((uint64_t)res * base) % mod);
    base = (uint32_t)(((uint64_t)base * base) % mod);
    exp >>= 1;
  }
  return res;
}

// fermat's little theorem
uint32_t mod_inverse(uint32_t a, uint32_t q) { 
  return mod_exp(a, q - 2, q); 
}

std::vector<uint32_t> generate_sequential_twiddles(uint32_t N, uint32_t q, uint32_t root) {
  std::vector<uint32_t> twiddles(N);
  for (uint32_t i = 0; i < N; i++) twiddles[i] = mod_exp(root, i, q);
  return twiddles;
}


/* ----------------------
    naive NTT O(N^2) 
------------------------- */

vector<uint32_t> naive_ntt(vector<uint32_t> a, uint32_t q, const vector<uint32_t>& W) {
    uint32_t N = a.size();
    vector<uint32_t> result(N, 0);

    for (uint32_t i = 0; i < N; i++) {
        uint32_t sum = 0;
        for (uint32_t j = 0; j < N; j++) {
            uint32_t idx = (i * j) % N;
            uint32_t coef = (uint32_t)(((uint64_t)a[j] * W[idx]) % q);
            sum = (uint32_t)((sum + coef) % q);
        }
        result[i] = sum;
    }
    return result;
}

vector<uint32_t> naive_intt(vector<uint32_t> a, uint32_t q, const vector<uint32_t> inv_W, const uint32_t inv_N) {
    a = naive_ntt(a, q, inv_W);
    for (uint32_t i = 0; i < a.size(); i++) {
      a[i] = (uint32_t)(((uint64_t)a[i] * inv_N) % q);
    }

    return a;
}

/*
------------------------
  fast NTT O(n log n)
------------------------
*/

std::vector<uint32_t> fast_gs_ntt(std::vector<uint32_t> a, uint32_t q, uint32_t root) {
  uint32_t N = a.size();
  for (uint32_t len = N; len >= 2; len >>= 1) {
    uint32_t wlen = mod_exp(root, N / len, q);
    for (uint32_t i = 0; i < N; i += len) {
      uint32_t w = 1;
      for (uint32_t j = 0; j < len / 2; j++) {
        uint32_t u = a[i + j];
        uint32_t v = a[i + j + len / 2];
        a[i + j] = (uint32_t)((u + v) % q);
        uint32_t diff = (uint32_t)((u + q - v) % q); 
        
        a[i + j + len / 2] = (uint32_t)(((uint64_t)diff * w) % q);
        w = (uint32_t)(((uint64_t)w * wlen) % q);
      }
    }
  }
  return a;
}

std::vector<uint32_t> fast_ct_intt(std::vector<uint32_t> a, uint32_t q, uint32_t inv_root, uint32_t inv_N) {
  uint32_t N = a.size();
  for (uint32_t len = 2; len <= N; len <<= 1) {
    uint32_t wlen = mod_exp(inv_root, N / len, q);
    for (uint32_t i = 0; i < N; i += len) {
      uint32_t w = 1;
      for (uint32_t j = 0; j < len / 2; j++) {
        uint32_t u = a[i + j];
        uint32_t v = (uint32_t)(((uint64_t)a[i + j + len / 2] * w) % q);
        
        a[i + j] = (uint32_t)((u + v) % q);
        a[i + j + len / 2] = (uint32_t)((u + q - v) % q);
        
        w = (uint32_t)(((uint64_t)w * wlen) % q);
      }
    }
  }
  for (uint32_t i = 0; i < N; i++) {
    a[i] = (uint32_t)(((uint64_t)a[i] * inv_N) % q);
  }
  return a;
}

/*
------------------------
  benchmaking
------------------------
*/

int main() {
  uint32_t N = 1 << 10;
  uint32_t cyclotomic_order = 2 * N;
  uint32_t bit_size = 30;
  uint32_t num_limbs = 16;

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
        // print first 10 coefficients of the first limb
        // if (i == 0 && j < 10){
        //     printf(" %lu ", original_rns_poly[i][j]);
        // }
    }
  }

  double naive_time = 0;
  double fast_time = 0;

  for (uint32_t i = 0; i < num_limbs; i++) {
    std::vector<uint32_t> poly = original_rns_poly[i];
    uint32_t q = rns_params[i].q;

    // generate twiddles
    std::vector<uint32_t> W = generate_sequential_twiddles(N, q, rns_params[i].root);
    std::vector<uint32_t> inv_W = generate_sequential_twiddles(N, q, rns_params[i].inv_root);

    // naive
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<uint32_t> res_naive = naive_ntt(poly, q, W);
    res_naive = naive_intt(res_naive, q, inv_W, rns_params[i].inv_N);
    auto end = std::chrono::high_resolution_clock::now();
    naive_time += std::chrono::duration<double, std::milli>(end - start).count();

    // fast
    start = std::chrono::high_resolution_clock::now();
    std::vector<uint32_t> res_fast = fast_gs_ntt(poly, q, rns_params[i].root);
    res_fast = fast_ct_intt(res_fast, q, rns_params[i].inv_root, rns_params[i].inv_N);
    end = std::chrono::high_resolution_clock::now();
    fast_time += std::chrono::duration<double, std::milli>(end - start).count();

    // correctness check
    for (uint32_t j = 0; j < N; j++) {
      if (res_naive[j] != original_rns_poly[i][j]) {
        printf("Naive NTT Error, mismatch at limb %d idx %d \n", i, j);
        printf("NTT Result: %u - Original Coefficient: %u \n", res_naive[j], original_rns_poly[i][j]);
        return 1;
      }
      if (res_fast[j] != original_rns_poly[i][j]) {
        printf("Fast NTT Error, mismatch at limb %d idx %d \n", i, j);
        printf("NTT Result: %u - Original Coefficient: %u \n", res_naive[j], original_rns_poly[i][j]);
        return 1;
      }
    }
  }

  printf("Naive NTT Time: %f ms\n", naive_time);
  printf("Fast NTT Time: %f ms\n", fast_time);
  return 0;

}