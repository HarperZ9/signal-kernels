// Shannon, Rényi, Tsallis, block, min, spectral, permutation entropy
// =============================================================================
// algorithms/entropy.hpp -- Information-theoretic entropy measures
//
// All functions operate on std::span<const T> for zero-copy consumption.
// Probabilities are assumed to be non-negative and sum to 1; counts are
// converted internally.  Log base 2 throughout (bits).
//
// Namespace: signal_kernels::algorithms
// =============================================================================

#pragma once
#ifndef SIGNAL_KERNELS_ENTROPY_HPP
#define SIGNAL_KERNELS_ENTROPY_HPP

#include "algorithms/_fft.hpp"
#include "algorithms/_numeric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_kernels::algorithms {

// ---------------------------------------------------------------------------
// shannon -- H = -sum(p * log2(p))  [bits]
// ---------------------------------------------------------------------------

[[nodiscard]] inline double shannon(std::span<const double> probs) noexcept {
    double h = 0.0;
    for (double p : probs) {
        if (p > 0.0) h -= p * std::log2(p);
    }
    return h;
}

[[nodiscard]] inline double
shannon_from_counts(std::span<const uint64_t> counts) noexcept {
    uint64_t total = 0;
    for (uint64_t c : counts) total += c;
    if (total == 0) return 0.0;
    double h = 0.0;
    const double inv = 1.0 / static_cast<double>(total);
    for (uint64_t c : counts) {
        if (c > 0) {
            double p = static_cast<double>(c) * inv;
            h -= p * std::log2(p);
        }
    }
    return h;
}

[[nodiscard]] inline double
shannon_from_bytes(std::span<const uint8_t> bytes) noexcept {
    uint64_t counts[256] = {};
    for (uint8_t b : bytes) ++counts[b];
    uint64_t total = bytes.size();
    if (total == 0) return 0.0;
    double h = 0.0;
    const double inv = 1.0 / static_cast<double>(total);
    for (size_t i = 0; i < 256; ++i) {
        if (counts[i] > 0) {
            double p = static_cast<double>(counts[i]) * inv;
            h -= p * std::log2(p);
        }
    }
    return h;
}

// ---------------------------------------------------------------------------
// renyi -- H_alpha = (1/(1-alpha)) * log2(sum(p^alpha))
//   alpha → 1 approaches Shannon; alpha = 2 is collision entropy.
// ---------------------------------------------------------------------------

[[nodiscard]] inline double renyi(std::span<const double> probs,
                                   double alpha) noexcept {
    if (alpha == 1.0) return shannon(probs);
    double sum = 0.0;
    for (double p : probs)
        if (p > 0.0) sum += std::pow(p, alpha);
    if (sum <= 0.0) return 0.0;
    return std::log2(sum) / (1.0 - alpha);
}

// ---------------------------------------------------------------------------
// tsallis -- S_q = (1 - sum(p^q)) / (q - 1)
//   q → 1 recovers Shannon (in nats*ln2 sense, here we keep bits via log2).
// ---------------------------------------------------------------------------

[[nodiscard]] inline double tsallis(std::span<const double> probs,
                                     double q) noexcept {
    if (std::abs(q - 1.0) < 1e-12) return shannon(probs);
    double sum = 0.0;
    for (double p : probs)
        if (p > 0.0) sum += std::pow(p, q);
    return (1.0 - sum) / (q - 1.0);
}

// ---------------------------------------------------------------------------
// min_entropy -- H_inf = -log2(max_p)
// ---------------------------------------------------------------------------

[[nodiscard]] inline double
min_entropy(std::span<const double> probs) noexcept {
    double mx = 0.0;
    for (double p : probs) mx = std::max(mx, p);
    return mx <= 0.0 ? 0.0 : -std::log2(mx);
}

// ---------------------------------------------------------------------------
// block_entropy -- entropy of non-overlapping L-grams over a byte sequence
//
// Partitions seq into consecutive L-byte blocks, counts their occurrences,
// and returns the Shannon entropy of the resulting distribution.
// ---------------------------------------------------------------------------

[[nodiscard]] inline double
block_entropy(std::span<const uint8_t> seq, size_t block_size) {
    if (seq.empty() || block_size == 0) return 0.0;
    if (block_size > seq.size()) return 0.0;

    std::map<std::vector<uint8_t>, uint64_t> counts;
    const size_t n_blocks = seq.size() / block_size;
    for (size_t i = 0; i < n_blocks; ++i) {
        std::vector<uint8_t> block(
            seq.data() + i * block_size,
            seq.data() + (i + 1) * block_size);
        ++counts[block];
    }

    std::vector<uint64_t> c;
    c.reserve(counts.size());
    for (auto& [k, v] : counts) c.push_back(v);
    return shannon_from_counts(c);
}

// ---------------------------------------------------------------------------
// spectral_entropy -- entropy of the normalized power spectrum
//
// Uses radix-2 FFT; pads signal to next power of two ≤ 2^16.
// Returns 0 if signal is too short or all power is concentrated.
// sample_rate unused in the entropy computation itself (kept for API symmetry).
// ---------------------------------------------------------------------------

[[nodiscard]] inline double
spectral_entropy(std::span<const double> signal,
                 [[maybe_unused]] double sample_rate = 1.0) {
    if (signal.size() < 2) return 0.0;

    // Find next power of two, cap at 2^16
    size_t n = 1;
    while (n < signal.size() && n < 65536) n <<= 1;
    if (n < signal.size()) n = 65536; // clamp

    std::vector<double> padded(n, 0.0);
    size_t copy_len = std::min(signal.size(), n);
    for (size_t i = 0; i < copy_len; ++i) padded[i] = signal[i];

    auto power = detail::power_spectrum(padded);

    // Normalize to probability distribution
    double total = 0.0;
    for (double p : power) total += p;
    if (total <= 0.0) return 0.0;

    std::vector<double> probs(power.size());
    for (size_t i = 0; i < power.size(); ++i)
        probs[i] = power[i] / total;

    return shannon(probs);
}

// ---------------------------------------------------------------------------
// permutation_entropy -- Bandt & Pompe (2002) ordinal pattern entropy
//
// order: embedding dimension (default 3), delay: time delay (default 1).
// Returns value in [0, log2(order!)] bits.
// ---------------------------------------------------------------------------

[[nodiscard]] inline double
permutation_entropy(std::span<const double> series,
                     int order = 3,
                     int delay = 1) noexcept {
    if (order < 2 || delay < 1) return 0.0;
    const size_t m   = static_cast<size_t>(order);
    const size_t tau = static_cast<size_t>(delay);
    const size_t n   = series.size();
    if (n < m * tau) return 0.0;

    // Encode each ordinal pattern as a rank-index tuple mapped to an integer
    // via factorial number system (Lehmer code).
    std::map<std::vector<uint8_t>, uint64_t> counts;
    std::vector<double>  embed(m);
    std::vector<uint8_t> pattern(m);

    const size_t limit = n - (m - 1) * tau;
    for (size_t i = 0; i < limit; ++i) {
        for (size_t j = 0; j < m; ++j)
            embed[j] = series[i + j * tau];

        // Compute rank permutation
        std::vector<size_t> idx(m);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
                  [&](size_t a, size_t b){ return embed[a] < embed[b]; });
        for (size_t j = 0; j < m; ++j)
            pattern[j] = static_cast<uint8_t>(idx[j]);

        ++counts[pattern];
    }

    std::vector<uint64_t> c;
    c.reserve(counts.size());
    for (auto& [k, v] : counts) c.push_back(v);
    return shannon_from_counts(c);
}

} // namespace signal_kernels::algorithms

#endif // SIGNAL_KERNELS_ENTROPY_HPP
