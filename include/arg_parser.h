#pragma once

#include <cstdint>
#include <string>

struct BenchConfig {
    uint32_t N = 1 << 10;
    uint32_t num_limbs = 32;
    uint32_t threads = 64;

    bool run_cpu_naive = false;
    bool run_cpu_fast = false;
    bool run_cpu_prod = false;
    bool run_gpu_radix2 = false;
    bool run_gpu_radix4 = false;

    bool use_barrett = false;

    // helper to check all test booleans
    bool run_all() const {
        return !run_cpu_naive && !run_cpu_fast && !run_cpu_prod && 
               !run_gpu_radix2 && !run_gpu_radix4;
    }
};

BenchConfig parse_args(int argc, char** argv);
