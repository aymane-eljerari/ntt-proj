#include "bench_utils.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>

bool correctness_check(const std::vector<std::vector<uint32_t>>& original, 
                       const std::vector<std::vector<uint32_t>>& computed) { 
    for (size_t i = 0; i < original.size(); i++) {
        if (original[i].size() != computed[i].size()) {
            return false;
        }
        
        for (size_t j = 0; j < original[i].size(); j++) {
            if (original[i][j] != computed[i][j]) {
                return false;
            }
        }
    }
    return true;
}

void run_benchmarks(const BenchConfig& config, 
                    const std::vector<RNSLimbParams>& rns_params,
                    const std::vector<std::vector<uint32_t>>& original_poly) {
    
    uint32_t N = 1 << config.N;
    uint32_t L = config.num_limbs;
    int num_runs = 5; 
    bool run_all = config.run_all();

    std::vector<BenchResult> results;

    // OpenMP warmup to avoid benchmark cold start
    int max_threads = omp_get_max_threads();
    std::vector<std::vector<uint32_t>> thread_warmup(max_threads);
    
    for (int i = 0; i < max_threads; i++) {
        thread_warmup[i].reserve(N); 
    }
    
    #pragma omp parallel for
    for (int i = 0; i < max_threads; i++) {
        thread_warmup[i].assign(N, 0); 
    }
    thread_warmup.clear(); 


    // GPU warmup
    if (run_all || config.run_gpu_radix2 || config.run_gpu_radix4) {
        std::vector<uint32_t> warmup_poly = original_poly[0];
        ntt_gpu_dif(warmup_poly, rns_params[0].q, rns_params[0].mu, rns_params[0].root, config.use_barrett);
    }

    // --- CPU Naive ---
    if (run_all || config.run_cpu_naive) {
        double total_time = 0;
        std::vector<std::vector<uint32_t>> final_result(L, std::vector<uint32_t>(N));
        
        for (int r = 0; r < num_runs; r++) {
            auto start = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < L; i++) {
                std::vector<uint32_t> poly = original_poly[i];
                std::vector<uint32_t> W = generate_sequential_twiddles(N, rns_params[i].q, rns_params[i].root);
                std::vector<uint32_t> inv_W = generate_sequential_twiddles(N, rns_params[i].q, rns_params[i].inv_root);
                
                std::vector<uint32_t> res;
                if (config.use_barrett) {
                    res = naive_ntt<true>(poly, rns_params[i].q, rns_params[i].mu, W);
                    res = naive_intt<true>(res, rns_params[i].q, rns_params[i].mu, inv_W, rns_params[i].inv_N);
                } else {
                    res = naive_ntt<false>(poly, rns_params[i].q, rns_params[i].mu, W);
                    res = naive_intt<false>(res, rns_params[i].q, rns_params[i].mu, inv_W, rns_params[i].inv_N);
                }
                
                if (r == num_runs - 1) final_result[i] = res;
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double, std::milli>(end - start).count();
        }
        bool passed = correctness_check(original_poly, final_result);
        results.push_back({"CPU Naive", total_time, total_time / num_runs, passed});
    }

    // --- CPU Fast ---
    if (run_all || config.run_cpu_fast) {
        double total_time = 0;
        std::vector<std::vector<uint32_t>> final_result(L, std::vector<uint32_t>(N));
        
        for (int r = 0; r < num_runs; r++) {
            auto start = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < L; i++) {
                std::vector<uint32_t> poly = original_poly[i];
                std::vector<uint32_t> res;
                
                if (config.use_barrett) {
                    res = fast_gs_ntt<true>(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].root);
                    res = fast_ct_intt<true>(res, rns_params[i].q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N);
                } else {
                    res = fast_gs_ntt<false>(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].root);
                    res = fast_ct_intt<false>(res, rns_params[i].q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N);
                }
                
                if (r == num_runs - 1) final_result[i] = res;
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double, std::milli>(end - start).count();
        }
        bool passed = correctness_check(original_poly, final_result);
        results.push_back({"CPU Fast", total_time, total_time / num_runs, passed});
    }

    // --- CPU Production ---
    if (run_all || config.run_cpu_prod) {
        double total_time = 0;
        std::vector<std::vector<uint32_t>> final_result(L, std::vector<uint32_t>(N));
        
        for (int r = 0; r < num_runs; r++) {
            auto start = std::chrono::high_resolution_clock::now();

            #pragma omp parallel for default(none) shared(L, original_poly, config, rns_params, final_result, r, num_runs)
            for (uint32_t i = 0; i < L; i++) {
                std::vector<uint32_t> poly = original_poly[i];
                std::vector<uint32_t> res;
                
                if (config.use_barrett) {
                    res = prod_gs_ntt<true>(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].omega_pow);
                    res = prod_ct_intt<true>(res, rns_params[i].q, rns_params[i].mu, rns_params[i].inv_omega_pow, rns_params[i].inv_N);
                } else {
                    res = prod_gs_ntt<false>(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].omega_pow);
                    res = prod_ct_intt<false>(res, rns_params[i].q, rns_params[i].mu, rns_params[i].inv_omega_pow, rns_params[i].inv_N);
                }
                
                if (r == num_runs - 1) final_result[i] = res;
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double, std::milli>(end - start).count();
        }
        bool passed = correctness_check(original_poly, final_result);
        results.push_back({"CPU Production", total_time, total_time / num_runs, passed});
    }

    // --- GPU Radix-2 ---
    if (run_all || config.run_gpu_radix2) {
        double total_time = 0;
        std::vector<std::vector<uint32_t>> final_result(L, std::vector<uint32_t>(N));
        
        for (int r = 0; r < num_runs; r++) {
            float run_time = 0;
            for (uint32_t i = 0; i < L; i++) {
                std::vector<uint32_t> poly = original_poly[i];
                run_time += ntt_gpu_dif(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].root, config.use_barrett);
                run_time += intt_gpu_dit(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N, config.use_barrett);
                
                if (r == num_runs - 1) final_result[i] = poly;
            }
            total_time += run_time;
        }
        bool passed = correctness_check(original_poly, final_result);
        results.push_back({"GPU Radix-2", total_time, total_time / num_runs, passed});
    }

    // --- GPU Radix-4 ---
    if (run_all || config.run_gpu_radix4) {
        double total_time = 0;
        std::vector<std::vector<uint32_t>> final_result(L, std::vector<uint32_t>(N));
        
        for (int r = 0; r < num_runs; r++) {
            float run_time = 0;
            for (uint32_t i = 0; i < L; i++) {
                std::vector<uint32_t> poly = original_poly[i];
                run_time += ntt_gpu_radix4_dif(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].root, config.use_barrett);
                run_time += intt_gpu_radix4_dit(poly, rns_params[i].q, rns_params[i].mu, rns_params[i].inv_root, rns_params[i].inv_N, config.use_barrett);
                
                if (r == num_runs - 1) final_result[i] = poly;
            }
            total_time += run_time;
        }
        bool passed = correctness_check(original_poly, final_result);
        results.push_back({"GPU Radix-4", total_time, total_time / num_runs, passed});
    }

    // --- OpenFHE OpenMP ---
    if (run_all || config.run_openfhe) {
        // get rid of lazy loading overhead
        int current_threads = omp_get_max_threads();
        omp_set_num_threads(1);
        benchmark_openfhe_rns_ntt(original_poly, N, 2 * N, rns_params, 1);
        
        // run the benchmark
        omp_set_num_threads(current_threads);
        OpenFHEBenchResult res = benchmark_openfhe_rns_ntt(original_poly, N, 2 * N, rns_params, num_runs);
        
        bool passed = correctness_check(original_poly, res.ntt_result);
        results.push_back({"OpenFHE (OMP)", res.time_ms * num_runs, res.time_ms, passed});
    }

    std::cout << "\n================ BENCHMARK RESULTS (Average over " << num_runs << " runs) ================\n";
    std::cout << std::left << std::setw(20) << "Implementation" 
              << " | " << std::setw(17) << "Avg Time (ms)" 
              << " | " << "Status\n";
    std::cout << "-----------------------------------------------------------------\n";
    for (const auto& res : results) {
        std::cout << std::left << std::setw(20) << res.name 
                  << " | " << std::fixed << std::setprecision(4) << std::setw(13) << res.avg_ms << " ms"
                  << " | " << (res.passed ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "=================================================================\n";

    if (!config.csv_out.empty()) {
        std::ifstream file_check(config.csv_out);
        bool write_header = !file_check.good();
        file_check.close();

        std::ofstream csv_file(config.csv_out, std::ios::app);
        if (csv_file.is_open()) {
            if (write_header) {
                csv_file << "Threads,N,L,Implementation,Total_Time_ms,Avg_Time_ms,Passed\n";
            }
            for (const auto& res : results) {
                csv_file << config.threads << ","
                         << N << ","
                         << L << ","
                         << res.name << ","
                         << res.total_ms << ","
                         << res.avg_ms << ","
                         << (res.passed ? "1" : "0") << "\n";
            }
            csv_file.close();
        } else {
            std::cerr << "Error: Could not open CSV file " << config.csv_out << " for appending.\n";
        }
    }
}

