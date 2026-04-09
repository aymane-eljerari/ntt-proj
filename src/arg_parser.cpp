#include "arg_parser.h"
#include <iostream>
#include <cstdlib>

BenchConfig parse_args(int argc, char** argv) {
    BenchConfig config;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-N" && i + 1 < argc) {
            config.N = std::stoi(argv[++i]);
        } else if (arg == "-L" && i + 1 < argc) {
            config.num_limbs = std::stoi(argv[++i]);
        } else if (arg == "-T" && i + 1 < argc) {
            config.threads = std::stoi(argv[++i]);
        } else if (arg == "--cpu-naive") {
            config.run_cpu_naive = true;
        } else if (arg == "--cpu-fast") {
            config.run_cpu_fast = true;
        } else if (arg == "--cpu-prod") {
            config.run_cpu_prod = true;
        } else if (arg == "--gpu-radix2") {
            config.run_gpu_radix2 = true;
        } else if (arg == "--gpu-radix4") {
            config.run_gpu_radix4 = true;
        } else if (arg == "--openfhe") {
            config.run_openfhe = true;
        } else if (arg == "--barrett") {
            config.use_barrett = true;
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: bench_ntt [options]\n"
                      << "Options:\n"
                      << "  -N <int>           Polynomial degree (default: 1024)\n"
                      << "  -L <int>           Number of limbs (default: 24)\n"
                      << "  -T <int>           Number of OpenFHE OpenMP threads (default: 1)\n"
                      << "  --cpu-naive        Run Naive CPU NTT\n"
                      << "  --cpu-fast         Run Fast CPU NTT\n"
                      << "  --cpu-prod         Run Production CPU NTT\n"
                      << "  --gpu-radix2       Run GPU Radix-2 NTT\n"
                      << "  --gpu-radix4       Run GPU Radix-4 NTT\n"
                      << "  --barrett          Enable Barrett reduction on GPU\n";
            std::exit(0);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            std::exit(1);
        }
    }

    return config;
}