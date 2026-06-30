// PELT change-point detection
//   Reference: Killick, Fearnhead, Eckley (2012) JASA
// =============================================================================
// algorithms/changepoint.hpp -- PELT change-point detection
//
// pelt() implements exact PELT (Pruned Exact Linear Time) segmentation.
// Cost functions: cost_l2, cost_l1, cost_poisson.
// Penalty helpers: bic_penalty, aic_penalty.
//
// Namespace: signal_kernels::algorithms
// =============================================================================

#pragma once
#ifndef SIGNAL_KERNELS_CHANGEPOINT_HPP
#define SIGNAL_KERNELS_CHANGEPOINT_HPP

#include "algorithms/_numeric.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <vector>

namespace signal_kernels::algorithms {

// ---------------------------------------------------------------------------
// Changepoint -- detected change point with location and segment cost
// ---------------------------------------------------------------------------

struct Changepoint {
    size_t index{0}; ///< First index of the new segment
    double cost{0.0};
};

// ---------------------------------------------------------------------------
// CostFn -- cost(series, start, end) over half-open interval [start, end)
//   Must satisfy subadditivity for PELT pruning to be exact.
// ---------------------------------------------------------------------------

using CostFn = std::function<double(std::span<const double>, size_t, size_t)>;

// ---------------------------------------------------------------------------
// cost_l2 -- sum of squared deviations from segment mean (L2 / variance cost)
//   C(s,t) = sum_{i=s}^{t-1} (x_i - mu)^2   where mu = mean(x[s..t))
//   Uses the identity: C = sum(x^2) - n*mu^2  via prefix sums.
// ---------------------------------------------------------------------------

namespace detail_cp {

struct PrefixSums {
    std::vector<double> sum;    // prefix sum of x
    std::vector<double> sum2;   // prefix sum of x^2
    explicit PrefixSums(std::span<const double> s)
        : sum(s.size() + 1, 0.0), sum2(s.size() + 1, 0.0) {
        for (size_t i = 0; i < s.size(); ++i) {
            sum[i + 1]  = sum[i]  + s[i];
            sum2[i + 1] = sum2[i] + s[i] * s[i];
        }
    }
    [[nodiscard]] double seg_sum(size_t lo, size_t hi) const noexcept {
        return sum[hi] - sum[lo];
    }
    [[nodiscard]] double seg_sum2(size_t lo, size_t hi) const noexcept {
        return sum2[hi] - sum2[lo];
    }
};

} // namespace detail_cp

[[nodiscard]] inline CostFn cost_l2_make(std::span<const double> series) {
    auto ps = std::make_shared<detail_cp::PrefixSums>(series);
    return [ps](std::span<const double> /*s*/, size_t lo, size_t hi) -> double {
        if (hi <= lo + 1) return 0.0;
        double n    = static_cast<double>(hi - lo);
        double s    = ps->seg_sum(lo, hi);
        double s2   = ps->seg_sum2(lo, hi);
        return s2 - (s * s) / n;
    };
}

// ---------------------------------------------------------------------------
// cost_l1 -- mean absolute deviation from median (L1 / robust cost)
//   Exact computation O((hi-lo) log(hi-lo)) per segment -- appropriate for
//   the penalty budget used in anti-cheat timing analysis.
// ---------------------------------------------------------------------------

[[nodiscard]] inline CostFn cost_l1_make(std::span<const double> series) {
    // Copy series for capture
    auto data = std::make_shared<std::vector<double>>(series.begin(),
                                                       series.end());
    return [data](std::span<const double> /*s*/, size_t lo, size_t hi) -> double {
        if (hi <= lo) return 0.0;
        std::vector<double> seg(data->begin() + static_cast<ptrdiff_t>(lo),
                                data->begin() + static_cast<ptrdiff_t>(hi));
        std::sort(seg.begin(), seg.end());
        double med = seg[seg.size() / 2];
        double cost = 0.0;
        for (double v : seg) cost += std::abs(v - med);
        return cost;
    };
}

// ---------------------------------------------------------------------------
// cost_poisson -- negative log-likelihood for Poisson rate model
//   C(s,t) = -2 * n * (mu * log(mu) - mu)  where mu = mean(x[s..t))
//   (drops constant terms that cancel in differences)
// ---------------------------------------------------------------------------

[[nodiscard]] inline CostFn cost_poisson_make(std::span<const double> series) {
    auto ps = std::make_shared<detail_cp::PrefixSums>(series);
    return [ps](std::span<const double> /*s*/, size_t lo, size_t hi) -> double {
        if (hi <= lo) return 0.0;
        double n   = static_cast<double>(hi - lo);
        double mu  = ps->seg_sum(lo, hi) / n;
        if (mu <= 0.0) return 0.0;
        return -2.0 * n * (mu * std::log(mu) - mu);
    };
}

// Global cost_l2 / cost_l1 / cost_poisson exposed as callable tag types
// so callers can pass them without constructing via make_* above.
// Usage: pelt(series, cost_l2, penalty, min_size)
// These are thin adapters that rebuild prefix sums per call -- use the
// _make() variants for repeated use on the same series.

namespace {
const CostFn cost_l2 = [](std::span<const double> s, size_t lo, size_t hi) {
    detail_cp::PrefixSums ps(s);
    if (hi <= lo + 1) return 0.0;
    double n  = static_cast<double>(hi - lo);
    double sv = ps.seg_sum(lo, hi);
    double s2 = ps.seg_sum2(lo, hi);
    return s2 - sv * sv / n;
};
const CostFn cost_l1 = [](std::span<const double> s, size_t lo, size_t hi) {
    if (hi <= lo) return 0.0;
    std::vector<double> seg(s.begin() + static_cast<ptrdiff_t>(lo),
                            s.begin() + static_cast<ptrdiff_t>(hi));
    std::sort(seg.begin(), seg.end());
    double med  = seg[seg.size() / 2];
    double cost = 0.0;
    for (double v : seg) cost += std::abs(v - med);
    return cost;
};
const CostFn cost_poisson = [](std::span<const double> s, size_t lo, size_t hi) {
    if (hi <= lo) return 0.0;
    double n = static_cast<double>(hi - lo);
    double mu = 0.0;
    for (size_t i = lo; i < hi; ++i) mu += s[i];
    mu /= n;
    if (mu <= 0.0) return 0.0;
    return -2.0 * n * (mu * std::log(mu) - mu);
};
} // anonymous namespace

// ---------------------------------------------------------------------------
// Penalty functions
// ---------------------------------------------------------------------------

[[nodiscard]] inline double bic_penalty(size_t n) noexcept {
    return std::log(static_cast<double>(n));
}

[[nodiscard]] inline double aic_penalty([[maybe_unused]] size_t n) noexcept {
    return 2.0;
}

// ---------------------------------------------------------------------------
// pelt -- Pruned Exact Linear Time change-point detection
//
// Killick, Fearnhead, Eckley (2012).
// Returns the list of detected change points (NOT including index 0 or n).
// Each Changepoint.index is the start of the new segment.
// min_size: minimum number of observations per segment (≥ 2).
// ---------------------------------------------------------------------------

[[nodiscard]] inline std::vector<Changepoint>
pelt(std::span<const double> series,
     const CostFn&           cost_fn,
     double                  penalty,
     size_t                  min_size = 2) {
    const size_t n = series.size();
    if (n < 2 * min_size) return {};

    // F[t] = optimal cost of segmenting series[0..t)
    // cp[t] = last change point before t in optimal segmentation
    std::vector<double> F(n + 1, std::numeric_limits<double>::infinity());
    std::vector<size_t> cp(n + 1, 0);
    F[0] = -penalty;

    std::vector<size_t> candidates = {0};

    for (size_t t = min_size; t <= n; ++t) {
        std::vector<size_t> next_candidates;
        double best = std::numeric_limits<double>::infinity();
        size_t best_s = 0;

        for (size_t s : candidates) {
            if (t - s < min_size) continue;
            double c = F[s] + cost_fn(series, s, t) + penalty;
            if (c < best) {
                best   = c;
                best_s = s;
            }
        }

        F[t]  = best;
        cp[t] = best_s;

        // Pruning: remove candidates that can never be optimal
        for (size_t s : candidates) {
            if (t - s < min_size) {
                next_candidates.push_back(s);
                continue;
            }
            double c_prune = F[s] + cost_fn(series, s, t);
            if (c_prune + penalty <= F[t]) {
                next_candidates.push_back(s);
            }
        }
        if (t >= min_size) next_candidates.push_back(t - min_size + 1);
        candidates = std::move(next_candidates);
    }

    // Backtrack
    std::vector<Changepoint> result;
    size_t cur = n;
    while (cp[cur] != 0) {
        size_t prev = cp[cur];
        result.push_back({prev, cost_fn(series, prev, cur)});
        cur = prev;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

} // namespace signal_kernels::algorithms

#endif // SIGNAL_KERNELS_CHANGEPOINT_HPP
