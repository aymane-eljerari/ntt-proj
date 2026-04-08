#include "ntt_openfhe.h"
#include <omp.h>
#include <chrono>

OpenFHEBenchResult benchmark_openfhe_rns_ntt(
    const std::vector<std::vector<uint32_t>>& original_rns_poly, 
    uint32_t N, 
    uint32_t cyclotomic_order,
    const std::vector<RNSLimbParams>& rns_params, 
    int num_runs) {

    uint32_t num_limbs = rns_params.size();

    std::vector<std::shared_ptr<ILNativeParams>> native_params;
    for (uint32_t i = 0; i < num_limbs; i++) {
        native_params.push_back(std::make_shared<ILNativeParams>(cyclotomic_order, rns_params[i].q, rns_params[i].root));
    }
    auto dcrt_params = std::make_shared<ILDCRTParams<BigInteger>>(cyclotomic_order, native_params);

    DCRTPoly openfhe_poly(dcrt_params, Format::COEFFICIENT);
    for (uint32_t i = 0; i < num_limbs; i++) {
        NativeVector vec(N, rns_params[i].q);
        for (uint32_t j = 0; j < N; j++) {
            vec[j] = original_rns_poly[i][j];
        }
        NativePoly np(native_params[i]);
        np.SetValues(vec, Format::COEFFICIENT);
        openfhe_poly.SetElementAtIndex(i, np);
    }

    // set and check OMP threads
    #ifdef BENCH_NTT
    int max_threads = omp_get_max_threads();
    int active_threads = 0;
    
    omp_set_num_threads(max_threads);
    
    #pragma omp parallel
    {
        #pragma omp single
        active_threads = omp_get_num_threads();
    }
    
    printf("\n--- OpenMP State ---\n");
    printf("Hardware Max Threads: %d\n", max_threads);
    printf("Threads spawned for testing: %d\n", active_threads);
    printf("Effective Thread Utilization: %d\n", std::min(active_threads, (int)num_limbs));
    printf("--------------------------\n\n");
    #endif

    double total_time = 0;
    DCRTPoly temp_poly;

    for (int run = 0; run < num_runs; run++) {
        temp_poly = openfhe_poly;

        auto start = std::chrono::high_resolution_clock::now();
        
        temp_poly.SetFormat(Format::EVALUATION);
        
        temp_poly.SetFormat(Format::COEFFICIENT);
        
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();
    }
 
    std::vector<std::vector<uint32_t>> extracted_result(num_limbs, std::vector<uint32_t>(N));
    for (uint32_t i = 0; i < num_limbs; i++) {
        NativePoly limb = openfhe_poly.GetElementAtIndex(i);
        for (uint32_t j = 0; j < N; j++) {
            extracted_result[i][j] = (uint32_t)limb[j].ConvertToInt();
        }
    }

    return { total_time / num_runs, extracted_result };
}