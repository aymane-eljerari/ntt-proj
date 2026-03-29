#pragma once

#include <vector>
#include <cstdint>

std::vector<uint32_t> naive_ntt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& W);
std::vector<uint32_t> naive_intt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t> inv_W, const uint32_t inv_N);

std::vector<uint32_t> fast_gs_ntt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, uint32_t root);
std::vector<uint32_t> fast_ct_intt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, uint32_t inv_root, uint32_t inv_N);

std::vector<uint32_t> prod_gs_ntt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& omega_pow);
std::vector<uint32_t> prod_ct_intt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& inv_omega_pow, uint32_t inv_N);