#pragma once

#include <vector>
#include <string>
#include <chrono>

#include "arg_parser.h"
#include "ntt_cpu.h"
#include "ntt_gpu.cuh"
#include "ntt_openfhe.h"

struct BenchResult {
    std::string name;
    double total_ms;
    double avg_ms;
    bool passed;
};

void run_benchmarks(const BenchConfig& config, 
                    const std::vector<RNSLimbParams>& rns_params,
                    const std::vector<std::vector<uint32_t>>& original_poly);


inline bool correctness_check(const std::vector<std::vector<uint32_t>>& original, 
                            const std::vector<std::vector<uint32_t>>& computed);