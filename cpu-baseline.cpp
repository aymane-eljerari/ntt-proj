#include <openfhe.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>

using namespace lbcrypto;
using namespace std;

/* ----------------------
    helper functions 
------------------------- */

uint64_t mod_exp(uint64_t base, uint64_t exp, uint64_t mod) {
  uint64_t res = 1;
  base = base % mod;
  while (exp > 0) {
    if (exp % 2 == 1) res = (uint64_t)((res * base) % mod);
    base = (uint64_t)((base * base) % mod);
    exp >>= 1;
  }
  return res;
}

// fermat's little theorem
uint64_t mod_inverse(uint64_t a, uint64_t q) { 
  return mod_exp(a, q - 2, q); 
}

struct RNSLimbParams {
  uint64_t q;
  uint64_t inv_N;
  uint64_t root;
  uint64_t inv_root;
};

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
  std::mt19937 gen(42);
  std::vector<std::vector<uint64_t>> original_rns_poly(num_limbs, std::vector<uint64_t>(N));
  for (uint32_t i = 0; i < num_limbs; i++) {
    std::uniform_int_distribution<uint64_t> dist(0, rns_params[i].q - 1);
    for (uint64_t j = 0; j < N; j++) {
        original_rns_poly[i][j] = dist(gen);
        // print first 10 coefficients of the first limb
        if (i == 0 && j < 10){
            printf(" %lu ", original_rns_poly[i][j]);
        }
    }
  }

  return 0;

}