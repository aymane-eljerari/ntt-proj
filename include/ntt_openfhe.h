#pragma once

#include <vector>
#include <cstdint>
#include <openfhe.h>
#include "utils.h"

using namespace lbcrypto;
using namespace std;

struct OpenFHEBenchResult {
    double time_ms;
    vector<vector<uint32_t>> ntt_result;
};


OpenFHEBenchResult benchmark_openfhe_rns_ntt(
    const vector<vector<uint32_t>>& original_rns_polynomial,
    uint32_t N,
    uint32_t cyclotomic_order,
    const vector<RNSLimbParams>& rns_params,
    int num_runs
); 
