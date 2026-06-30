// Granger causality
// =============================================================================
// algorithms/causal.hpp -- Granger causality test
//
// granger_causality(x, y, max_lag):
//   Tests whether x "Granger-causes" y.  Fits a restricted AR(lag) on y and
//   an unrestricted ARX(lag) on y with x lags, then computes the F-statistic.
//   Selects the lag with the best (largest) F statistic up to max_lag.
//
// GrangerResult: { f_stat, p_value, optimal_lag }
//
// p-value approximation: Snedecor F CDF via regularized incomplete Beta.
//
// Namespace: signal_kernels::algorithms
// =============================================================================

#pragma once
#ifndef SIGNAL_KERNELS_CAUSAL_HPP
#define SIGNAL_KERNELS_CAUSAL_HPP

#include "algorithms/_numeric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace signal_kernels::algorithms {

// ---------------------------------------------------------------------------
// GrangerResult
// ---------------------------------------------------------------------------

struct GrangerResult {
    double f_stat{0.0};
    double p_value{1.0};
    int    optimal_lag{0};
};

namespace detail_causal {

// Regularized incomplete beta I_x(a, b) via continued fraction (Lentz method)
// For computing F-distribution CDF: p = I_{df2/(df2+df1*x)}(df2/2, df1/2)
[[nodiscard]] inline double
reg_incomplete_beta(double x, double a, double b) noexcept {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    // Use symmetry if x > (a+1)/(a+b+2)
    bool sym = (x > (a + 1.0) / (a + b + 2.0));
    if (sym) return 1.0 - reg_incomplete_beta(1.0 - x, b, a);

    // Log prefactor
    double lbeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
    double front = std::exp(std::log(x) * a + std::log(1.0 - x) * b - lbeta) / a;

    // Modified Lentz continued fraction
    double f = 1.0, C = 1.0, D = 0.0;
    for (int m = 0; m <= 200; ++m) {
        for (int sign = 0; sign <= 1; ++sign) {
            double d;
            if (sign == 0 && m == 0) {
                d = 1.0;
            } else if (sign == 0) {
                d = -(a + m - 1.0) * (a + b + m - 1.0) * x
                    / ((a + 2.0 * m - 2.0) * (a + 2.0 * m - 1.0));
            } else {
                d = static_cast<double>(m) * (b - m) * x
                    / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
            }
            D = 1.0 + d * D;
            if (std::abs(D) < 1e-30) D = 1e-30;
            D = 1.0 / D;
            C = 1.0 + d / C;
            if (std::abs(C) < 1e-30) C = 1e-30;
            f *= C * D;
            if (std::abs(C * D - 1.0) < 1e-8) break;
        }
    }
    return std::clamp(front * (f - 1.0), 0.0, 1.0);
}

[[nodiscard]] inline double f_pvalue(double F, double df1, double df2) noexcept {
    if (F <= 0.0 || df1 <= 0.0 || df2 <= 0.0) return 1.0;
    double x = df2 / (df2 + df1 * F);
    // P(F_stat > F) = I_x(df2/2, df1/2)
    return reg_incomplete_beta(x, df2 / 2.0, df1 / 2.0);
}

// OLS residual sum of squares for y ~ X (no intercept, caller prepends column of ones)
[[nodiscard]] inline double
ols_ssr(const std::vector<std::vector<double>>& X,
        const std::vector<double>& y) {
    // X is [n x p]; solve normal equations X'X * beta = X'y
    const size_t n = y.size();
    const size_t p = X.empty() ? 0 : X[0].size();
    if (n == 0 || p == 0) return 0.0;

    SmallMatrix A(p, p, 0.0);
    std::vector<double> b(p, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t r = 0; r < p; ++r) {
            b[r] += X[i][r] * y[i];
            for (size_t c = 0; c < p; ++c)
                A.at(r, c) += X[i][r] * X[i][c];
        }
    }
    if (!A.solve_in_place(b)) return 0.0; // singular
    // beta is now in b; compute residuals
    double ssr = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double pred = 0.0;
        for (size_t r = 0; r < p; ++r) pred += X[i][r] * b[r];
        double e = y[i] - pred;
        ssr += e * e;
    }
    return ssr;
}

} // namespace detail_causal

// ---------------------------------------------------------------------------
// granger_causality -- tests if x Granger-causes y
// ---------------------------------------------------------------------------

[[nodiscard]] inline GrangerResult
granger_causality(std::span<const double> x,
                   std::span<const double> y,
                   int max_lag) noexcept {
    if (max_lag < 1) return {};
    const size_t n = std::min(x.size(), y.size());
    if (n <= static_cast<size_t>(2 * max_lag + 2)) return {};

    GrangerResult best;
    best.f_stat = -1.0;

    for (int lag = 1; lag <= max_lag; ++lag) {
        const size_t eff = n - static_cast<size_t>(lag);
        if (eff < 3) continue;

        // Build design matrices
        // Restricted: y_t = c + sum_{j=1}^{lag} a_j * y_{t-j}
        // Unrestricted: same + sum_{j=1}^{lag} b_j * x_{t-j}
        std::vector<std::vector<double>> Xr(eff, std::vector<double>(static_cast<size_t>(lag) + 1));
        std::vector<std::vector<double>> Xu(eff, std::vector<double>(2u * static_cast<size_t>(lag) + 1));
        std::vector<double> yt(eff);

        for (size_t i = 0; i < eff; ++i) {
            size_t t = i + static_cast<size_t>(lag);
            yt[i] = y[t];
            Xr[i][0] = 1.0; // intercept
            Xu[i][0] = 1.0;
            for (int j = 0; j < lag; ++j) {
                Xr[i][static_cast<size_t>(j) + 1] = y[t - static_cast<size_t>(j) - 1];
                Xu[i][static_cast<size_t>(j) + 1] = y[t - static_cast<size_t>(j) - 1];
                Xu[i][static_cast<size_t>(lag + j) + 1] = x[t - static_cast<size_t>(j) - 1];
            }
        }

        double ssr_r = detail_causal::ols_ssr(Xr, yt);
        double ssr_u = detail_causal::ols_ssr(Xu, yt);

        double df1 = static_cast<double>(lag);
        double df2 = static_cast<double>(eff) - 2.0 * static_cast<double>(lag) - 1.0;

        if (df2 <= 0.0 || ssr_u <= 0.0) continue;

        double F = ((ssr_r - ssr_u) / df1) / (ssr_u / df2);
        if (F < 0.0) F = 0.0;

        if (F > best.f_stat) {
            best.f_stat     = F;
            best.optimal_lag = lag;
            best.p_value    = detail_causal::f_pvalue(F, df1, df2);
        }
    }
    if (best.f_stat < 0.0) best.f_stat = 0.0;
    return best;
}

} // namespace signal_kernels::algorithms

#endif // SIGNAL_KERNELS_CAUSAL_HPP
