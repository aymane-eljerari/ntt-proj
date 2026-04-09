#pragma once
#include <vector>
#include <cstdint>

// naive
template <bool UseBarrett>
std::vector<uint32_t> naive_ntt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& W);
template <bool UseBarrett>
std::vector<uint32_t> naive_intt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& inv_W, uint32_t inv_N);

// fast
template <bool UseBarrett>
std::vector<uint32_t> fast_gs_ntt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, uint32_t root);
template <bool UseBarrett>
std::vector<uint32_t> fast_ct_intt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, uint32_t inv_root, uint32_t inv_N);

// prod
template <bool UseBarrett>
std::vector<uint32_t> prod_gs_ntt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& omega_pow);
template <bool UseBarrett>
std::vector<uint32_t> prod_ct_intt(std::vector<uint32_t> a, uint32_t q, uint64_t mu, const std::vector<uint32_t>& inv_omega_pow, uint32_t inv_N);
