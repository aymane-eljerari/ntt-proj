#include <openfhe.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <omp.h>

#include "utils.h"
#include "ntt_cpu.h"
#include "ntt_openfhe.h"
#include "ntt_gpu.cuh"
#include "arg_parser.h"
#include "bench_utils.h"

using namespace lbcrypto;
using namespace std;

int main(int argc, char** argv) {
    BenchConfig config = parse_args(argc, argv);
    
    uint32_t N = 1 << config.N;
    uint32_t num_limbs = config.num_limbs;
    uint32_t cyclotomic_order = 2 * N;
    uint32_t bit_size = 28;

    omp_set_num_threads(config.threads);

    printf("================ CONFIGURATION ================\n");
    printf("Degree (N)      : %u\n", N);
    printf("Limbs (L)       : %u\n", num_limbs);
    printf("OMP Threads     : %u\n", config.threads);
    printf("Barrett GPU     : %s\n", config.use_barrett ? "ON" : "OFF");
    printf("===============================================\n\n");

    // OpenFHE RNS Polynomial setup
    std::vector<RNSLimbParams> rns_params(num_limbs);
    NativeInteger current_prime = FirstPrime<NativeInteger>(bit_size, cyclotomic_order);

    for (uint32_t i = 0; i < num_limbs; i++) {
        NativeInteger psi_unity = RootOfUnity<NativeInteger>(cyclotomic_order, current_prime);
        uint32_t q = (uint32_t)current_prime.ConvertToInt();
        
        uint32_t omega = mod_exp<false>((uint32_t)psi_unity.ConvertToInt(), 2, q);
        
        uint32_t inv_omega = mod_inverse(omega, q);
        
        uint32_t inv_N = mod_inverse(N, q);

        rns_params[i].q = q;

        // precomputed barett reduction factor
        rns_params[i].mu = (uint64_t)(((__uint128_t) 1 << 64) / q);
        rns_params[i].root = omega;
        rns_params[i].inv_root = inv_omega;
        rns_params[i].inv_N = inv_N;

        // twiddle precomputations
        rns_params[i].omega_pow = generate_sequential_twiddles(N, q, omega);
        rns_params[i].inv_omega_pow = generate_sequential_twiddles(N, q, inv_omega);
        
        current_prime = NextPrime<NativeInteger>(current_prime, cyclotomic_order);
    }

    // populate RNS polynomials
    mt19937 gen(42);
    std::vector<std::vector<uint32_t>> original_rns_poly(num_limbs, vector<uint32_t>(N));
    for (uint32_t i = 0; i < num_limbs; i++) {
        uniform_int_distribution<uint32_t> dist(0, rns_params[i].q - 1);
        for (uint32_t j = 0; j < N; j++) {
            original_rns_poly[i][j] = dist(gen);
        }
    }

    run_benchmarks(config, rns_params, original_rns_poly);

    return 0;
}