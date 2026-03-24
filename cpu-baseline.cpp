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
  uint64_t q;
  uint64_t inv_N;
  uint64_t root;
  uint64_t inv_root;
};


/* ----------------------
    helper functions 
------------------------- */

uint64_t mod_exp(uint64_t base, uint64_t exp, uint64_t mod) {
  uint64_t res = 1;
  base = base % mod;
  while (exp > 0) {
    if (exp % 2 == 1) res = (uint64_t)(((unsigned __int128)res * base) % mod);
    base = (uint64_t)(((unsigned __int128)base * base) % mod);
    exp >>= 1;
  }
  return res;
}

// fermat's little theorem
uint64_t mod_inverse(uint64_t a, uint64_t q) { 
  return mod_exp(a, q - 2, q); 
}

vector<uint64_t> generate_sequential_twiddles(uint64_t N, uint64_t q, uint64_t root) {
  vector<uint64_t> twiddles(N);
  for (uint64_t i = 0; i < N; i++) twiddles[i] = mod_exp(root, i, q);
  return twiddles;
}


/* ----------------------
    naive NTT O(N^2) 
------------------------- */

vector<uint64_t> naive_ntt(vector<uint64_t> a, uint64_t q, const vector<uint64_t>& W) {
    uint64_t N = a.size();
    vector<uint64_t> result(N, 0);

    for (uint64_t i = 0; i < N; i++) {
        uint64_t sum = 0;
        for (uint64_t j = 0; j < N; j++) {
            uint64_t idx = (i * j) % N;
            uint64_t coef = (uint64_t)(((unsigned __int128)a[j] * W[idx]) % q);
            sum = (sum + coef) % q;
        }
        result[i] = sum;
    }
    return result;
}

vector<uint64_t> naive_intt(vector<uint64_t> a, uint64_t q, const vector<uint64_t> inv_W, const uint64_t inv_N) {
    a = naive_ntt(a, q, inv_W);
    for (uint64_t i = 0; i < a.size(); i++) {
      a[i] = (uint64_t)(((unsigned __int128)a[i] * inv_N) % q);
    }

    return a;
}

int main() {
  uint32_t N = 1 << 12;
  uint32_t cyclotomic_order = 2 * N;
  uint32_t bit_size = 54;
  uint32_t num_limbs = 32;

  // init rns polynomial + get the first prime
  std::vector<RNSLimbParams> rns_params(num_limbs);
  NativeInteger current_prime = FirstPrime<NativeInteger>(bit_size, cyclotomic_order);

  // generate and store rns parameters
  for (uint32_t i = 0; i < num_limbs; i++) {
    NativeInteger psi_unity = RootOfUnity<NativeInteger>(cyclotomic_order, current_prime);
    uint64_t q = current_prime.ConvertToInt();
    uint64_t omega = mod_exp(psi_unity.ConvertToInt(), 2, q);
    uint64_t inv_omega = mod_inverse(omega, q);
    uint64_t inv_N = mod_inverse(N, q);

    rns_params[i].q = q;
    rns_params[i].root = omega;
    rns_params[i].inv_root = inv_omega;
    rns_params[i].inv_N = inv_N;

    current_prime = NextPrime<NativeInteger>(current_prime, cyclotomic_order);

  }

  // rng sampling to populate the rns limbs
  mt19937 gen(42);
  std::vector<std::vector<uint64_t>> original_rns_poly(num_limbs, vector<uint64_t>(N));
  for (uint32_t i = 0; i < num_limbs; i++) {
    uniform_int_distribution<uint64_t> dist(0, rns_params[i].q - 1);
    for (uint64_t j = 0; j < N; j++) {
        original_rns_poly[i][j] = dist(gen);
        // print first 10 coefficients of the first limb
        // if (i == 0 && j < 10){
        //     printf(" %lu ", original_rns_poly[i][j]);
        // }
    }
  }

  double naive_time = 0;

  for (uint32_t i = 0; i < num_limbs; i++) {
    std::vector<uint64_t> poly = original_rns_poly[i];
    uint64_t q = rns_params[i].q;

    // generate twiddle
    std::vector<uint64_t> W = generate_sequential_twiddles(N, q, rns_params[i].root);
    std::vector<uint64_t> inv_W = generate_sequential_twiddles(N, q, rns_params[i].inv_root);

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<uint64_t> res_naive = naive_ntt(poly, q, W);
    res_naive = naive_intt(res_naive, q, inv_W, rns_params[i].inv_N);
    auto end = std::chrono::high_resolution_clock::now();
    naive_time += std::chrono::duration<double, std::milli>(end - start).count();

    // correctness check
    for (uint32_t j = 0; j < N; j++) {
      if (res_naive[j] != original_rns_poly[i][j]) {
        printf("Error, mismatch at limb %d idx %d \n", i, j);
        printf("NTT Result: %ld - Original Coefficient: %ld \n", res_naive[j], original_rns_poly[i][j]);
        return 1;
      }
    }
  }

  

  printf("Naive NTT completed successfully\n");
  printf("Naive NTT Time: %f ms\n", naive_time);
  return 0;

}