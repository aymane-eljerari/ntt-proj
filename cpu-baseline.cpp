#include <openfhe.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>

using namespace lbcrytpo;


/* helper functions */
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

int main() {
  uint32_t N = 1 << 12;
  uint32_t cyclotomic_order = 2 * N;
  uint32_t bit_size = 50;
  uint32_t num_limbs = 5;

}