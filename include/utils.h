#pragma once

#include <cstdint>
#include <vector>

struct RNSLimbParams {
    uint32_t q;
    uint64_t mu;
    uint32_t inv_N;
    uint32_t root;
    uint32_t inv_root;
    std::vector<uint32_t> omega_pow;
    std::vector<uint32_t> inv_omega_pow;
};

#ifdef USE_BARRETT

inline uint32_t mod_add(uint32_t a, uint32_t b, uint32_t q) {
    uint32_t sum = a + b;
    return (sum >= q) ? sum - q : sum;
}

inline uint32_t mod_sub(uint32_t a, uint32_t b, uint32_t q) {
    uint32_t diff = a + q - b;
    return (diff >= q) ? diff - q : diff;
}

inline uint32_t mod_mul(uint32_t a, uint32_t b, uint32_t q, uint64_t mu) {
    uint64_t ab = (uint64_t)a * b;
    uint64_t q1 = ((__uint128_t)ab * mu) >> 64;
    uint64_t r = ab - q1 * q;
    return (r >= q) ? (uint32_t)(r - q) : (uint32_t)r;
}

#else

inline uint32_t mod_add(uint32_t a, uint32_t b, uint32_t q) { return (a + b) % q; }
inline uint32_t mod_sub(uint32_t a, uint32_t b, uint32_t q) { return (a + q - b) % q; }
inline uint32_t mod_mul(uint32_t a, uint32_t b, uint32_t q, uint64_t mu) { return ((uint64_t)a * b) % q; }

#endif

inline uint32_t mod_exp(uint32_t base, uint32_t exp, uint32_t mod) {
    uint32_t res = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (uint32_t)(((uint64_t)res * base) % mod);
        base = (uint32_t)(((uint64_t)base * base) % mod);
        exp >>= 1;
    }
    return res;
}

inline uint32_t mod_inverse(uint32_t a, uint32_t q) { 
    return mod_exp(a, q - 2, q); 
}

inline std::vector<uint32_t> generate_sequential_twiddles(uint32_t N, uint32_t q, uint32_t root) {
    std::vector<uint32_t> twiddles(N);
    for (uint32_t i = 0; i < N; i++) twiddles[i] = mod_exp(root, i, q);
    return twiddles;
}