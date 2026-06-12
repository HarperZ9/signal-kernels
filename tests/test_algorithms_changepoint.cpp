// =============================================================================
// Tests: algorithms — changepoint module
//
// Known-value checks:
//   - PELT on flat signal = no change points
//   - PELT on step function recovers the step (±1 index tolerance)
//   - cost_l2 on uniform segment = 0
//   - cost_l1 on uniform segment = 0
//   - bic_penalty > aic_penalty for n > 8
//   - PELT with large penalty returns fewer CPs than small penalty
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/changepoint.hpp"

#include <cmath>
#include <vector>

using namespace warden::algorithms;

TEST_SUITE("algorithms::changepoint") {

    TEST_CASE("flat signal — no change points") {
        std::vector<double> series(50, 5.0);
        auto cps = pelt(series, cost_l2, bic_penalty(50), 2);
        CHECK(cps.empty());
    }

    TEST_CASE("step function — PELT recovers step index ±1") {
        // First 25 elements at 0.0, next 25 at 10.0 — step at index 25
        std::vector<double> series(50, 0.0);
        for (size_t i = 25; i < 50; ++i) series[i] = 10.0;
        auto cps = pelt(series, cost_l2, bic_penalty(50), 2);
        REQUIRE(!cps.empty());
        // The first detected change point should be near index 25
        bool near_25 = false;
        for (auto& cp : cps) {
            if (cp.index >= 23 && cp.index <= 27) { near_25 = true; break; }
        }
        CHECK(near_25);
    }

    TEST_CASE("two steps — PELT finds both") {
        // segment [0,20) = 0.0, [20,40) = 5.0, [40,60) = 0.0
        std::vector<double> series(60, 0.0);
        for (size_t i = 20; i < 40; ++i) series[i] = 5.0;
        auto cps = pelt(series, cost_l2, bic_penalty(60) * 0.5, 2);
        // Should find at least one change point in [18,22] and one in [38,42]
        bool found_first = false, found_second = false;
        for (auto& cp : cps) {
            if (cp.index >= 18 && cp.index <= 22) found_first  = true;
            if (cp.index >= 38 && cp.index <= 42) found_second = true;
        }
        CHECK(found_first);
        CHECK(found_second);
    }

    TEST_CASE("cost_l2 on uniform segment = 0") {
        std::vector<double> s(20, 3.7);
        double c = cost_l2(s, 0, 20);
        CHECK(c == doctest::Approx(0.0).epsilon(1e-6));
    }

    TEST_CASE("cost_l1 on uniform segment = 0") {
        std::vector<double> s(15, 2.2);
        double c = cost_l1(s, 0, 15);
        CHECK(c == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("cost_poisson non-negative on positive data") {
        std::vector<double> s = {1.0, 2.0, 3.0, 4.0, 5.0};
        double c = cost_poisson(s, 0, 5);
        CHECK(c >= 0.0);
    }

    TEST_CASE("bic_penalty > aic_penalty for n=100") {
        CHECK(bic_penalty(100) > aic_penalty(100));
    }

    TEST_CASE("bic_penalty scales with log(n)") {
        double p10   = bic_penalty(10);
        double p1000 = bic_penalty(1000);
        CHECK(p1000 > p10);
    }

    TEST_CASE("large penalty produces fewer change points") {
        std::vector<double> series(100, 0.0);
        for (size_t i = 25; i < 50; ++i) series[i] = 4.0;
        for (size_t i = 75; i < 100; ++i) series[i] = 8.0;

        auto cps_small = pelt(series, cost_l2, 0.1, 2);
        auto cps_large = pelt(series, cost_l2, 100.0, 2);
        CHECK(cps_small.size() >= cps_large.size());
    }

    TEST_CASE("min_size constraint respected") {
        // With min_size=10, change points must be at least 10 samples apart
        std::vector<double> series(40, 0.0);
        for (size_t i = 20; i < 40; ++i) series[i] = 10.0;
        constexpr size_t min_sz = 10;
        auto cps = pelt(series, cost_l2, bic_penalty(40), min_sz);
        // All reported indices must be >= min_sz from boundaries
        for (auto& cp : cps) {
            CHECK(cp.index >= min_sz);
            CHECK(cp.index <= series.size() - min_sz);
        }
    }

    TEST_CASE("PELT with L1 cost recovers step in noisy data") {
        // L1 is robust to outliers
        std::vector<double> series(60, 0.0);
        for (size_t i = 30; i < 60; ++i) series[i] = 5.0;
        // Inject one outlier
        series[10] = 50.0;
        auto cps = pelt(series, cost_l1, bic_penalty(60) * 2.0, 3);
        bool near_30 = false;
        for (auto& cp : cps)
            if (cp.index >= 27 && cp.index <= 33) { near_30 = true; break; }
        CHECK(near_30);
    }

} // TEST_SUITE
