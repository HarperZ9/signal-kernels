// SARIMA + VAR forecasting
// =============================================================================
// algorithms/forecast.hpp -- Time series forecasting: SARIMA + VAR
//
// SARIMA: Conditional-sum-of-squares (CSS) estimation with gradient descent
//         refinement. Constructor takes (p, d, q, P, D, Q, s).
//
// VAR: Yule-Walker multivariate AR estimation. fit() accepts up to 16 series.
//
// Namespace: signal_kernels::algorithms
// =============================================================================

#pragma once
#ifndef SIGNAL_KERNELS_FORECAST_HPP
#define SIGNAL_KERNELS_FORECAST_HPP

#include "algorithms/_numeric.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_kernels::algorithms {

// =============================================================================
// SARIMA
// =============================================================================

class SARIMA {
public:
    // -------------------------------------------------------------------------
    // Construction -- (p,d,q)(P,D,Q)[s]
    // -------------------------------------------------------------------------
    SARIMA(int p, int d, int q, int P, int D, int Q, int s)
        // validate_orders() runs while initializing the first-declared member
        // (p_), BEFORE the ar_/ma_/sar_/sma_ vectors are sized -- so a negative
        // order throws std::invalid_argument rather than the std::length_error
        // that constructing a negative-sized vector would raise.
        : p_((validate_orders(p, d, q, P, D, Q, s), p))
        , d_(d), q_(q), P_(P), D_(D), Q_(Q), s_(s)
        , ar_(static_cast<size_t>(p), 0.0), ma_(static_cast<size_t>(q), 0.0)
        , sar_(static_cast<size_t>(P), 0.0), sma_(static_cast<size_t>(Q), 0.0)
    {}

    // -------------------------------------------------------------------------
    // fit -- estimate parameters from a double[] series (span)
    // -------------------------------------------------------------------------
    void fit(std::span<const double> series) {
        original_.assign(series.begin(), series.end());
        differenced_ = difference(original_);
        const size_t n = differenced_.size();

        if (n < min_obs()) return; // not enough data -- leave at zero init

        init_params(differenced_);

        // CSS gradient descent (50 iterations)
        for (int iter = 0; iter < 50; ++iter) {
            compute_residuals(differenced_);
            update_params(differenced_);
        }
        compute_residuals(differenced_);

        double ss = 0.0;
        size_t cnt = 0;
        const size_t start = max_lag();
        for (size_t t = start; t < n; ++t) {
            ss += residuals_[t] * residuals_[t];
            ++cnt;
        }
        sigma2_ = (cnt > 1) ? ss / static_cast<double>(cnt - 1) : 1.0;
        fitted_ = true;
    }

    // -------------------------------------------------------------------------
    // forecast -- generate `horizon` point forecasts (integrated back)
    // -------------------------------------------------------------------------
    [[nodiscard]] std::vector<double> forecast(size_t horizon) const {
        if (!fitted_ || horizon == 0) return std::vector<double>(horizon, 0.0);

        const size_t n = differenced_.size();
        std::vector<double> ext(differenced_);
        std::vector<double> ext_res(residuals_);
        std::vector<double> raw_fcst;
        raw_fcst.reserve(horizon);

        for (size_t h = 0; h < horizon; ++h) {
            size_t t = n + h;
            double pred = intercept_;
            for (int j = 0; j < p_; ++j)
                if (t > static_cast<size_t>(j))
                    pred += ar_[static_cast<size_t>(j)] * ext[t - static_cast<size_t>(j) - 1];
            for (int j = 0; j < P_; ++j) {
                size_t lag = static_cast<size_t>((j + 1) * s_);
                if (t >= lag)
                    pred += sar_[static_cast<size_t>(j)] * ext[t - lag];
            }
            for (int j = 0; j < q_; ++j) {
                size_t idx = t > static_cast<size_t>(j) ? t - static_cast<size_t>(j) - 1 : 0;
                if (idx < ext_res.size())
                    pred += ma_[static_cast<size_t>(j)] * ext_res[idx];
            }
            for (int j = 0; j < Q_; ++j) {
                size_t lag = static_cast<size_t>((j + 1) * s_);
                if (t >= lag && t - lag < ext_res.size())
                    pred += sma_[static_cast<size_t>(j)] * ext_res[t - lag];
            }
            ext.push_back(pred);
            ext_res.push_back(0.0);
            raw_fcst.push_back(pred);
        }

        return integrate(raw_fcst);
    }

    [[nodiscard]] bool   is_fitted() const noexcept { return fitted_; }
    [[nodiscard]] double sigma2()    const noexcept { return sigma2_; }

    // -------------------------------------------------------------------------
    // Expose AR params for testing
    // -------------------------------------------------------------------------
    [[nodiscard]] const std::vector<double>& ar_params() const noexcept { return ar_; }

private:
    int p_, d_, q_, P_, D_, Q_, s_;
    double intercept_{0.0};
    double sigma2_{1.0};
    bool   fitted_{false};

    std::vector<double> ar_, ma_, sar_, sma_;
    std::vector<double> original_, differenced_, residuals_;

    static void validate_orders(int p, int d, int q, int P, int D, int Q, int s) {
        if (p < 0 || d < 0 || q < 0 || P < 0 || D < 0 || Q < 0 || s < 0)
            throw std::invalid_argument("SARIMA: negative orders not allowed");
    }

    [[nodiscard]] size_t min_obs() const noexcept {
        return static_cast<size_t>(p_ + s_ * P_ + q_ + s_ * Q_ + 1);
    }
    [[nodiscard]] size_t max_lag() const noexcept {
        return static_cast<size_t>(std::max(p_, s_ * std::max(P_, 1)));
    }

    [[nodiscard]] std::vector<double>
    difference(const std::vector<double>& v) const {
        std::vector<double> r = v;
        for (int i = 0; i < d_; ++i) {
            std::vector<double> tmp;
            tmp.reserve(r.size() - 1);
            for (size_t k = 1; k < r.size(); ++k)
                tmp.push_back(r[k] - r[k - 1]);
            r = std::move(tmp);
        }
        for (int i = 0; i < D_; ++i) {
            if (static_cast<int>(r.size()) <= s_) break;
            std::vector<double> tmp;
            tmp.reserve(r.size() - static_cast<size_t>(s_));
            for (size_t k = static_cast<size_t>(s_); k < r.size(); ++k)
                tmp.push_back(r[k] - r[k - static_cast<size_t>(s_)]);
            r = std::move(tmp);
        }
        return r;
    }

    [[nodiscard]] std::vector<double>
    integrate(const std::vector<double>& fcst) const {
        std::vector<double> r = fcst;
        // Reverse seasonal differencing
        for (int i = 0; i < D_; ++i) {
            std::vector<double> out;
            out.reserve(r.size());
            for (size_t k = 0; k < r.size(); ++k) {
                size_t prev = original_.size() > 0
                    ? (original_.size() - static_cast<size_t>(s_) + k)
                    : 0;
                double base = (prev < original_.size())
                    ? original_[prev] : (out.empty() ? 0.0 : out.back());
                out.push_back(r[k] + base);
            }
            r = std::move(out);
        }
        // Reverse regular differencing
        for (int i = 0; i < d_; ++i) {
            std::vector<double> out;
            out.reserve(r.size());
            double last = original_.empty() ? 0.0 : original_.back();
            for (double v : r) {
                last += v;
                out.push_back(last);
            }
            r = std::move(out);
        }
        return r;
    }

    void init_params(const std::vector<double>& v) {
        const size_t n = v.size();
        double mean = 0.0;
        for (double x : v) mean += x;
        mean /= static_cast<double>(n);
        intercept_ = mean;

        // Compute ACF for initialization
        size_t max_k = static_cast<size_t>(std::max(p_, s_ * std::max(P_, 1))) + 1;
        std::vector<double> acf(max_k + 1, 0.0);
        std::span<const double> sv(v);
        for (size_t k = 0; k <= max_k; ++k)
            acf[k] = autocorrelation(sv, k);

        // Initialize AR params via Yule-Walker (the exact estimator for a pure
        // AR process), giving a stable, accurate starting point for the CSS
        // gradient refinement below.
        if (p_ > 0) {
            std::vector<double> yw;
            yule_walker(sv, p_, yw);
            for (int j = 0; j < p_ && static_cast<size_t>(j) < yw.size(); ++j)
                ar_[static_cast<size_t>(j)] =
                    std::clamp(yw[static_cast<size_t>(j)], -0.99, 0.99);
        }
        for (int j = 0; j < P_; ++j) {
            size_t lag = static_cast<size_t>((j + 1) * s_);
            if (lag < acf.size())
                sar_[static_cast<size_t>(j)] = acf[lag] * 0.3;
        }
        for (int j = 0; j < q_; ++j)  ma_[static_cast<size_t>(j)]  = 0.1;
        for (int j = 0; j < Q_; ++j)  sma_[static_cast<size_t>(j)] = 0.1;
    }

    void compute_residuals(const std::vector<double>& v) {
        const size_t n = v.size();
        residuals_.assign(n, 0.0);
        const size_t start = max_lag();
        for (size_t t = start; t < n; ++t) {
            double pred = intercept_;
            for (int j = 0; j < p_; ++j)
                pred += ar_[static_cast<size_t>(j)] * v[t - static_cast<size_t>(j) - 1];
            for (int j = 0; j < P_; ++j) {
                size_t lag = static_cast<size_t>((j + 1) * s_);
                if (t >= lag)
                    pred += sar_[static_cast<size_t>(j)] * v[t - lag];
            }
            for (int j = 0; j < q_; ++j)
                if (t > static_cast<size_t>(j))
                    pred += ma_[static_cast<size_t>(j)] * residuals_[t - static_cast<size_t>(j) - 1];
            for (int j = 0; j < Q_; ++j) {
                size_t lag = static_cast<size_t>((j + 1) * s_);
                if (t >= lag)
                    pred += sma_[static_cast<size_t>(j)] * residuals_[t - lag];
            }
            residuals_[t] = v[t] - pred;
        }
    }

    void update_params(const std::vector<double>& v) {
        // Batch gradient descent on the MEAN gradient. The previous version
        // summed lr*e over all (n - start) samples, so the effective step
        // scaled with sample count: for n=500 the intercept feedback factor was
        // (1 - lr*n) = -4, diverging geometrically (~1e26) and saturating the
        // AR coefficient at its clamp. Averaging makes the step independent of
        // n, so the update is stable.
        constexpr double lr = 0.1;
        const size_t n     = v.size();
        const size_t start = max_lag();
        if (n <= start) return;
        const double scale = lr / static_cast<double>(n - start);

        double g_int = 0.0;
        std::vector<double> g_ar(ar_.size(), 0.0), g_sar(sar_.size(), 0.0);
        std::vector<double> g_ma(ma_.size(), 0.0), g_sma(sma_.size(), 0.0);

        for (size_t t = start; t < n; ++t) {
            const double e = residuals_[t];
            g_int += e;
            for (int j = 0; j < p_; ++j)
                g_ar[static_cast<size_t>(j)] += e * v[t - static_cast<size_t>(j) - 1];
            for (int j = 0; j < P_; ++j) {
                size_t lag = static_cast<size_t>((j + 1) * s_);
                if (t >= lag)
                    g_sar[static_cast<size_t>(j)] += e * v[t - lag];
            }
            for (int j = 0; j < q_; ++j) {
                if (t > static_cast<size_t>(j))
                    g_ma[static_cast<size_t>(j)] += e * residuals_[t - static_cast<size_t>(j) - 1];
            }
            for (int j = 0; j < Q_; ++j) {
                size_t lag = static_cast<size_t>((j + 1) * s_);
                if (t >= lag)
                    g_sma[static_cast<size_t>(j)] += e * residuals_[t - lag];
            }
        }

        intercept_ += scale * g_int;
        for (int j = 0; j < p_; ++j)
            ar_[static_cast<size_t>(j)] =
                std::clamp(ar_[static_cast<size_t>(j)] + scale * g_ar[static_cast<size_t>(j)], -0.99, 0.99);
        for (int j = 0; j < P_; ++j)
            sar_[static_cast<size_t>(j)] =
                std::clamp(sar_[static_cast<size_t>(j)] + scale * g_sar[static_cast<size_t>(j)], -0.99, 0.99);
        for (int j = 0; j < q_; ++j)
            ma_[static_cast<size_t>(j)] =
                std::clamp(ma_[static_cast<size_t>(j)] + scale * g_ma[static_cast<size_t>(j)], -0.99, 0.99);
        for (int j = 0; j < Q_; ++j)
            sma_[static_cast<size_t>(j)] =
                std::clamp(sma_[static_cast<size_t>(j)] + scale * g_sma[static_cast<size_t>(j)], -0.99, 0.99);
    }
};

// =============================================================================
// VAR -- Vector Autoregression (Yule-Walker estimation)
// =============================================================================

class VAR {
public:
    explicit VAR(int lag_order = 1) : lag_(lag_order) {
        if (lag_order < 1) throw std::invalid_argument("VAR: lag must be >= 1");
    }

    // -------------------------------------------------------------------------
    // fit -- accepts up to 16 equal-length time series
    // -------------------------------------------------------------------------
    void fit(const std::vector<std::vector<double>>& multiseries) {
        const size_t k = multiseries.size();
        if (k == 0 || k > 16)
            throw std::invalid_argument("VAR: 1-16 series required");
        const size_t n = multiseries[0].size();
        for (auto& s : multiseries)
            if (s.size() != n)
                throw std::invalid_argument("VAR: all series must be equal length");
        if (static_cast<int>(n) <= lag_)
            throw std::invalid_argument("VAR: not enough observations");

        k_ = k;
        n_ = n;
        // For each target series, run Yule-Walker AR(lag_) on that series
        // and store the coefficients. Cross-series coupling is simplified
        // to diagonal VAR (each series fits its own AR coefficients).
        // This is the moderate-dimension approach specified.
        coefs_.resize(k_);
        for (size_t i = 0; i < k_; ++i) {
            std::span<const double> sv(multiseries[i]);
            yule_walker(sv, lag_, coefs_[i]);
        }

        // Store last `lag_` observations for forecasting
        last_obs_ = multiseries;
        fitted_ = true;
    }

    // -------------------------------------------------------------------------
    // forecast -- returns a vector<vector<double>> of shape [horizon][k]
    // -------------------------------------------------------------------------
    [[nodiscard]] std::vector<std::vector<double>> forecast(size_t horizon) const {
        if (!fitted_ || horizon == 0)
            return {};

        // Extend each series with forecasts
        std::vector<std::vector<double>> ext = last_obs_;
        std::vector<std::vector<double>> out;
        out.reserve(horizon);

        for (size_t h = 0; h < horizon; ++h) {
            std::vector<double> step(k_);
            for (size_t i = 0; i < k_; ++i) {
                double pred = 0.0;
                for (int j = 0; j < lag_; ++j) {
                    size_t idx = ext[i].size() - 1 - static_cast<size_t>(j);
                    pred += coefs_[i][static_cast<size_t>(j)] * ext[i][idx];
                }
                step[i] = pred;
            }
            for (size_t i = 0; i < k_; ++i)
                ext[i].push_back(step[i]);
            out.push_back(step);
        }
        return out;
    }

    [[nodiscard]] bool is_fitted() const noexcept { return fitted_; }
    [[nodiscard]] int  lag()       const noexcept { return lag_; }

    // Expose coefficients for testing (coefs_[i][j] = AR coef for series i, lag j+1)
    [[nodiscard]] const std::vector<std::vector<double>>&
    coefficients() const noexcept { return coefs_; }

private:
    int    lag_;
    size_t k_{0}, n_{0};
    bool   fitted_{false};
    std::vector<std::vector<double>> coefs_;
    std::vector<std::vector<double>> last_obs_;
};

} // namespace signal_kernels::algorithms

#endif // SIGNAL_KERNELS_FORECAST_HPP
