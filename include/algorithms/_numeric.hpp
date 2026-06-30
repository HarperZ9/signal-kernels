// Numeric utilities
// =============================================================================
// algorithms/_numeric.hpp -- Internal numeric primitives
//
// Welford streaming variance, log-sum-exp, banded/small-matrix linear algebra.
// All routines are template-free or constrained to double; stdlib only.
//
// Namespace: signal_kernels::algorithms
// =============================================================================

#pragma once
#ifndef SIGNAL_KERNELS_NUMERIC_HPP
#define SIGNAL_KERNELS_NUMERIC_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_kernels::algorithms {

// ---------------------------------------------------------------------------
// WelfordAccumulator -- O(1)-per-sample numerically stable variance
// ---------------------------------------------------------------------------

class WelfordAccumulator {
public:
    void push(double x) noexcept {
        ++count_;
        double delta  = x - mean_;
        mean_        += delta / static_cast<double>(count_);
        double delta2 = x - mean_;
        m2_          += delta * delta2;
    }

    [[nodiscard]] double mean()     const noexcept { return mean_; }
    [[nodiscard]] double variance() const noexcept {
        return count_ < 2 ? 0.0 : m2_ / static_cast<double>(count_ - 1);
    }
    [[nodiscard]] double population_variance() const noexcept {
        return count_ < 1 ? 0.0 : m2_ / static_cast<double>(count_);
    }
    [[nodiscard]] double stddev()  const noexcept { return std::sqrt(variance()); }
    [[nodiscard]] uint64_t count() const noexcept { return count_; }
    void reset() noexcept { count_ = 0; mean_ = 0.0; m2_ = 0.0; }

private:
    uint64_t count_{0};
    double   mean_{0.0};
    double   m2_{0.0};
};

// ---------------------------------------------------------------------------
// log_sum_exp -- numerically stable log(sum(exp(x_i)))
// ---------------------------------------------------------------------------

[[nodiscard]] inline double log_sum_exp(std::span<const double> xs) noexcept {
    if (xs.empty()) return -std::numeric_limits<double>::infinity();
    double mx = *std::max_element(xs.begin(), xs.end());
    if (!std::isfinite(mx)) return mx;
    double sum = 0.0;
    for (double v : xs) sum += std::exp(v - mx);
    return mx + std::log(sum);
}

// ---------------------------------------------------------------------------
// SmallMatrix -- dense row-major matrix ≤ 32x32, heap-allocated via vector
//
// Used for Yule-Walker and SARIMA normal equations.
// ---------------------------------------------------------------------------

class SmallMatrix {
public:
    SmallMatrix() = default;
    SmallMatrix(size_t rows, size_t cols, double fill = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, fill) {}

    [[nodiscard]] double& at(size_t r, size_t c) noexcept {
        return data_[r * cols_ + c];
    }
    [[nodiscard]] const double& at(size_t r, size_t c) const noexcept {
        return data_[r * cols_ + c];
    }
    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }

    // In-place Gauss-Jordan to solve A*x = b (overwrites *this and b).
    // Returns false if singular.
    [[nodiscard]] bool solve_in_place(std::vector<double>& b) {
        const size_t n = rows_;
        if (n == 0 || n != cols_ || n != b.size()) return false;

        // Augment
        SmallMatrix aug(n, n + 1);
        for (size_t r = 0; r < n; ++r) {
            for (size_t c = 0; c < n; ++c)
                aug.at(r, c) = at(r, c);
            aug.at(r, n) = b[r];
        }

        for (size_t col = 0; col < n; ++col) {
            // Partial pivot
            size_t pivot = col;
            for (size_t r = col + 1; r < n; ++r)
                if (std::abs(aug.at(r, col)) > std::abs(aug.at(pivot, col)))
                    pivot = r;
            if (pivot != col)
                for (size_t c = 0; c <= n; ++c)
                    std::swap(aug.at(col, c), aug.at(pivot, c));

            double diag = aug.at(col, col);
            if (std::abs(diag) < 1e-15) return false;

            for (size_t c = col; c <= n; ++c)
                aug.at(col, c) /= diag;

            for (size_t r = 0; r < n; ++r) {
                if (r == col) continue;
                double factor = aug.at(r, col);
                for (size_t c = col; c <= n; ++c)
                    aug.at(r, c) -= factor * aug.at(col, c);
            }
        }

        for (size_t r = 0; r < n; ++r)
            b[r] = aug.at(r, n);
        return true;
    }

private:
    size_t rows_{0}, cols_{0};
    std::vector<double> data_;
};

// ---------------------------------------------------------------------------
// autocorrelation -- sample ACF at lag k (Pearson normalization)
// ---------------------------------------------------------------------------

[[nodiscard]] inline double autocorrelation(std::span<const double> series,
                                             size_t k) noexcept {
    const size_t n = series.size();
    if (n <= k + 1) return 0.0;
    double mean = 0.0;
    for (double v : series) mean += v;
    mean /= static_cast<double>(n);

    double denom = 0.0, numer = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = series[i] - mean;
        denom += d * d;
    }
    for (size_t i = k; i < n; ++i)
        numer += (series[i] - mean) * (series[i - k] - mean);

    return denom < 1e-15 ? 0.0 : numer / denom;
}

// ---------------------------------------------------------------------------
// build_toeplitz_rhs -- helpers for Yule-Walker
// Build the R*phi = r system: R is Toeplitz from acf[0..p-1], r is acf[1..p]
// ---------------------------------------------------------------------------

inline void yule_walker(std::span<const double> series,
                        int p,
                        std::vector<double>& phi) {
    if (p <= 0) { phi.clear(); return; }
    const size_t up = static_cast<size_t>(p);

    std::vector<double> acf(up + 1);
    for (size_t k = 0; k <= up; ++k)
        acf[k] = autocorrelation(series, k);

    // Toeplitz solve via Durbin-Levinson is O(p^2) but SmallMatrix is fine ≤32
    SmallMatrix R(up, up);
    std::vector<double> rhs(up);
    for (size_t i = 0; i < up; ++i) {
        rhs[i] = acf[i + 1];
        for (size_t j = 0; j < up; ++j) {
            size_t lag = (i >= j) ? (i - j) : (j - i);
            R.at(i, j) = acf[lag];
        }
    }

    phi.resize(up, 0.0);
    if (!R.solve_in_place(rhs)) {
        std::fill(phi.begin(), phi.end(), 0.0);
        return;
    }
    phi = rhs;
}

} // namespace signal_kernels::algorithms

#endif // SIGNAL_KERNELS_NUMERIC_HPP
