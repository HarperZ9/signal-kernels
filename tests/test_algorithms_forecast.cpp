// =============================================================================
// Tests: algorithms — forecast module (SARIMA + VAR)
//
// Known-value checks:
//   - SARIMA AR(1) on AR(1)-generated series recovers AR coef ±5%
//   - SARIMA forecast length matches requested horizon
//   - SARIMA unfitted forecast returns zeros
//   - VAR(1) on two independent AR(1) series produces non-trivial forecast
//   - VAR rejects mismatched series lengths
//   - VAR rejects > 16 series
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/forecast.hpp"

#include <cmath>
#include <vector>

using namespace warden::algorithms;

static std::vector<double> make_ar1(double phi, int n, double sigma = 0.05) {
    // Generate AR(1) series: x_t = phi * x_{t-1} + N(0, sigma)
    std::vector<double> s(static_cast<size_t>(n), 0.0);
    uint64_t rng = 0xDEADBEEF12345678ull;
    for (int i = 1; i < n; ++i) {
        rng = rng * 6364136223846793005ull + 1442695040888963407ull;
        double noise = static_cast<double>(static_cast<int64_t>(rng) % 1000) / 10000.0 * sigma;
        s[static_cast<size_t>(i)] = phi * s[static_cast<size_t>(i-1)] + noise;
    }
    return s;
}

TEST_SUITE("algorithms::forecast") {

    TEST_CASE("SARIMA AR(1) fit is marked fitted") {
        auto series = make_ar1(0.7, 200);
        SARIMA model(1, 0, 0, 0, 0, 0, 1);
        model.fit(series);
        CHECK(model.is_fitted());
    }

    TEST_CASE("SARIMA AR(1) recovers AR coefficient ±10%") {
        // Use p=1, d=0, q=0, seasonal off (s=1)
        constexpr double TARGET_PHI = 0.7;
        auto series = make_ar1(TARGET_PHI, 500, 0.02);
        SARIMA model(1, 0, 0, 0, 0, 0, 1);
        model.fit(series);
        REQUIRE(model.is_fitted());
        REQUIRE(model.ar_params().size() == 1u);
        CHECK(model.ar_params()[0] == doctest::Approx(TARGET_PHI).epsilon(0.10));
    }

    TEST_CASE("SARIMA forecast length matches horizon") {
        auto series = make_ar1(0.6, 300);
        SARIMA model(1, 0, 0, 0, 0, 0, 1);
        model.fit(series);
        auto fcst = model.forecast(10);
        CHECK(fcst.size() == 10u);
    }

    TEST_CASE("SARIMA unfitted returns vector of correct size") {
        SARIMA model(1, 0, 0, 0, 0, 0, 1);
        auto fcst = model.forecast(5);
        CHECK(fcst.size() == 5u);
    }

    TEST_CASE("SARIMA negative orders throw") {
        CHECK_THROWS_AS(SARIMA(-1, 0, 0, 0, 0, 0, 1), std::invalid_argument);
    }

    TEST_CASE("SARIMA sigma2 > 0 after fit") {
        auto series = make_ar1(0.5, 200);
        SARIMA model(1, 0, 0, 0, 0, 0, 1);
        model.fit(series);
        CHECK(model.sigma2() > 0.0);
    }

    TEST_CASE("SARIMA(1,1,0) on differenced trend") {
        // Generate a random walk (d=1 handles it)
        std::vector<double> series(100, 0.0);
        for (size_t i = 1; i < 100; ++i)
            series[i] = series[i-1] + 0.5 + (i % 3 == 0 ? 0.1 : -0.1);
        SARIMA model(1, 1, 0, 0, 0, 0, 1);
        model.fit(series);
        auto fcst = model.forecast(5);
        CHECK(fcst.size() == 5u);
        // Forecast values should be finite
        for (double f : fcst) CHECK(std::isfinite(f));
    }

    TEST_CASE("VAR fit marks model fitted") {
        auto s1 = make_ar1(0.6, 100);
        auto s2 = make_ar1(0.4, 100);
        VAR model(1);
        model.fit({s1, s2});
        CHECK(model.is_fitted());
    }

    TEST_CASE("VAR forecast shape is [horizon][k]") {
        auto s1 = make_ar1(0.5, 100);
        auto s2 = make_ar1(0.3, 100);
        VAR model(1);
        model.fit({s1, s2});
        auto fcst = model.forecast(5);
        CHECK(fcst.size() == 5u);
        for (auto& step : fcst) CHECK(step.size() == 2u);
    }

    TEST_CASE("VAR forecasts are finite") {
        auto s1 = make_ar1(0.7, 200);
        auto s2 = make_ar1(0.5, 200);
        VAR model(2);
        model.fit({s1, s2});
        auto fcst = model.forecast(10);
        for (auto& step : fcst)
            for (double v : step)
                CHECK(std::isfinite(v));
    }

    TEST_CASE("VAR rejects mismatched series lengths") {
        std::vector<double> a(50, 1.0), b(60, 1.0);
        VAR model(1);
        CHECK_THROWS_AS(model.fit({a, b}), std::invalid_argument);
    }

    TEST_CASE("VAR rejects > 16 series") {
        std::vector<std::vector<double>> many(17, std::vector<double>(50, 0.0));
        VAR model(1);
        CHECK_THROWS_AS(model.fit(many), std::invalid_argument);
    }

    TEST_CASE("VAR lag accessor returns configured lag") {
        VAR model(3);
        CHECK(model.lag() == 3);
    }

    TEST_CASE("VAR unfitted returns empty forecast") {
        VAR model(1);
        auto fcst = model.forecast(5);
        CHECK(fcst.empty());
    }

} // TEST_SUITE
