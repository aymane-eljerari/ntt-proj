#pragma once
#include <vector>
#include <cstdint>

void ntt_gpu_dif(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t root_of_unity);
void intt_gpu_dit(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t inv_root_of_unity, uint32_t inv_N);

void ntt_gpu_radix4_dif(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t root_of_unity);
void intt_gpu_radix4_dit(std::vector<uint32_t>& data, uint32_t q, uint64_t mu, uint32_t inv_root_of_unity, uint32_t inv_N);