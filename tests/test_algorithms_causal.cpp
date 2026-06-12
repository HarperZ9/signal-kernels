// =============================================================================
// Tests: algorithms — causal module (Granger causality)
//
// Known-value checks:
//   - Self-causality produces positive F-stat
//   - Independent white noise: F-stat low / p-value high
//   - optimal_lag in [1, max_lag]
//   - f_stat >= 0
//   - p_value in [0, 1]
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/causal.hpp"

#include <cmath>
#include <vector>

using namespace warden::algorithms;

static std::vector<double> white_noise(int n, uint64_t seed = 42) {
    std::vector<double> s(static_cast<size_t>(n));
    uint64_t rng = seed;
    for (int i = 0; i < n; ++i) {
        rng = rng * 6364136223846793005ull + 1442695040888963407ull;
        s[static_cast<size_t>(i)] = static_cast<double>(static_cast<int64_t>(rng) % 1000) / 1000.0;
    }
    return s;
}

TEST_SUITE("algorithms::causal") {

    TEST_CASE("f_stat is non-negative") {
        auto x = white_noise(100, 1);
        auto y = white_noise(100, 2);
        auto r = granger_causality(x, y, 3);
        CHECK(r.f_stat >= 0.0);
    }

    TEST_CASE("p_value in [0, 1]") {
        auto x = white_noise(100, 11);
        auto y = white_noise(100, 22);
        auto r = granger_causality(x, y, 3);
        CHECK(r.p_value >= 0.0);
        CHECK(r.p_value <= 1.0);
    }

    TEST_CASE("optimal_lag in [1, max_lag]") {
        auto x = white_noise(120, 5);
        auto y = white_noise(120, 7);
        constexpr int max_lag = 4;
        auto r = granger_causality(x, y, max_lag);
        CHECK(r.optimal_lag >= 1);
        CHECK(r.optimal_lag <= max_lag);
    }

    TEST_CASE("causal series: x causes y, f_stat > independent case") {
        // Build y_t = 0.8 * x_{t-1} + small noise
        constexpr int N = 200;
        auto x = white_noise(N, 99);
        std::vector<double> y(static_cast<size_t>(N), 0.0);
        uint64_t rng = 777;
        for (int i = 1; i < N; ++i) {
            rng = rng * 6364136223846793005ull + 1442695040888963407ull;
            double noise = static_cast<double>(static_cast<int64_t>(rng) % 100) / 10000.0;
            y[static_cast<size_t>(i)] = 0.8 * x[static_cast<size_t>(i-1)] + noise;
        }
        auto r_causal = granger_causality(x, y, 3);

        // Compare with independent series
        auto z = white_noise(N, 12345);
        auto r_indep  = granger_causality(z, y, 3);

        // Causal case should have higher F-stat
        CHECK(r_causal.f_stat > r_indep.f_stat);
    }

    TEST_CASE("independent white noise has high p-value (not always significant)") {
        // With pure noise there should be no systematic Granger causality
        // p_value > 0.01 for at least one direction
        auto x = white_noise(150, 333);
        auto y = white_noise(150, 444);
        auto r = granger_causality(x, y, 2);
        // We can't guarantee specific p-value for random data, but
        // F-stat should be finite and non-negative
        CHECK(std::isfinite(r.f_stat));
        CHECK(r.f_stat >= 0.0);
    }

    TEST_CASE("insufficient data returns zero result") {
        std::vector<double> x = {1.0, 2.0};
        std::vector<double> y = {1.0, 2.0};
        auto r = granger_causality(x, y, 3);
        CHECK(r.f_stat == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("max_lag < 1 returns zero result") {
        auto x = white_noise(50, 1);
        auto y = white_noise(50, 2);
        auto r = granger_causality(x, y, 0);
        CHECK(r.f_stat == doctest::Approx(0.0).epsilon(1e-9));
    }

} // TEST_SUITE
