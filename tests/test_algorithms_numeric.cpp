// =============================================================================
// Tests: algorithms — _numeric module
//
// Known-value checks:
//   - WelfordAccumulator matches naive variance
//   - log_sum_exp numerically stable
//   - SmallMatrix Gauss-Jordan solve
//   - autocorrelation lag-0 = 1.0
//   - yule_walker AR(1) recovers coefficient
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/_numeric.hpp"

#include <cmath>
#include <vector>

using namespace warden::algorithms;

TEST_SUITE("algorithms::_numeric") {

    TEST_CASE("WelfordAccumulator mean and variance") {
        WelfordAccumulator w;
        std::vector<double> vals = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
        for (double v : vals) w.push(v);
        CHECK(w.mean()     == doctest::Approx(5.0).epsilon(1e-9));
        CHECK(w.variance() == doctest::Approx(4.571428571).epsilon(1e-6));
        CHECK(w.count()    == 8u);
    }

    TEST_CASE("WelfordAccumulator stddev matches sqrt(variance)") {
        WelfordAccumulator w;
        for (int i = 1; i <= 10; ++i) w.push(static_cast<double>(i));
        CHECK(w.stddev() == doctest::Approx(std::sqrt(w.variance())).epsilon(1e-12));
    }

    TEST_CASE("WelfordAccumulator reset") {
        WelfordAccumulator w;
        w.push(1.0); w.push(2.0); w.push(3.0);
        w.reset();
        CHECK(w.count() == 0u);
        CHECK(w.mean()  == doctest::Approx(0.0).epsilon(1e-12));
    }

    TEST_CASE("WelfordAccumulator single element variance = 0") {
        WelfordAccumulator w;
        w.push(42.0);
        CHECK(w.variance() == doctest::Approx(0.0).epsilon(1e-12));
    }

    TEST_CASE("log_sum_exp empty returns -inf") {
        std::vector<double> empty;
        double r = log_sum_exp(empty);
        CHECK(std::isinf(r));
        CHECK(r < 0.0);
    }

    TEST_CASE("log_sum_exp single element returns that element") {
        std::vector<double> v = {3.7};
        CHECK(log_sum_exp(v) == doctest::Approx(3.7).epsilon(1e-12));
    }

    TEST_CASE("log_sum_exp numerically stable large values") {
        // log(exp(1000) + exp(1000)) = 1000 + log(2)
        std::vector<double> v = {1000.0, 1000.0};
        CHECK(log_sum_exp(v) == doctest::Approx(1000.0 + std::log(2.0)).epsilon(1e-9));
    }

    TEST_CASE("SmallMatrix Gauss-Jordan 2x2 identity") {
        SmallMatrix A(2, 2);
        A.at(0,0) = 1.0; A.at(0,1) = 0.0;
        A.at(1,0) = 0.0; A.at(1,1) = 1.0;
        std::vector<double> b = {3.0, 7.0};
        bool ok = A.solve_in_place(b);
        CHECK(ok);
        CHECK(b[0] == doctest::Approx(3.0).epsilon(1e-9));
        CHECK(b[1] == doctest::Approx(7.0).epsilon(1e-9));
    }

    TEST_CASE("SmallMatrix Gauss-Jordan 3x3 known") {
        // 2x + y - z = 8
        // -3x - y + 2z = -11
        // -2x + y + 2z = -3
        // Solution: x=2, y=3, z=-1
        SmallMatrix A(3, 3);
        A.at(0,0) =  2; A.at(0,1) =  1; A.at(0,2) = -1;
        A.at(1,0) = -3; A.at(1,1) = -1; A.at(1,2) =  2;
        A.at(2,0) = -2; A.at(2,1) =  1; A.at(2,2) =  2;
        std::vector<double> b = {8.0, -11.0, -3.0};
        bool ok = A.solve_in_place(b);
        CHECK(ok);
        CHECK(b[0] == doctest::Approx(2.0).epsilon(1e-9));
        CHECK(b[1] == doctest::Approx(3.0).epsilon(1e-9));
        CHECK(b[2] == doctest::Approx(-1.0).epsilon(1e-9));
    }

    TEST_CASE("SmallMatrix singular returns false") {
        SmallMatrix A(2, 2);
        A.at(0,0) = 1.0; A.at(0,1) = 2.0;
        A.at(1,0) = 2.0; A.at(1,1) = 4.0;
        std::vector<double> b = {1.0, 2.0};
        bool ok = A.solve_in_place(b);
        CHECK(!ok);
    }

    TEST_CASE("autocorrelation lag 0 = 1.0") {
        std::vector<double> s = {1.0, 2.0, 3.0, 4.0, 5.0};
        CHECK(autocorrelation(s, 0) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("autocorrelation lag beyond series returns 0") {
        std::vector<double> s = {1.0, 2.0};
        CHECK(autocorrelation(s, 5) == doctest::Approx(0.0).epsilon(1e-12));
    }

    TEST_CASE("yule_walker AR(1) recovers coefficient") {
        // Generate AR(1) series: x_t = 0.8 * x_{t-1} + epsilon
        constexpr int N = 1000;
        std::vector<double> series(N, 0.0);
        uint64_t rng = 123456789ull;
        for (int i = 1; i < N; ++i) {
            rng = rng * 6364136223846793005ull + 1442695040888963407ull;
            double noise = static_cast<double>(static_cast<int64_t>(rng) % 10000) / 100000.0;
            series[static_cast<size_t>(i)] = 0.8 * series[static_cast<size_t>(i-1)] + noise;
        }
        std::vector<double> phi;
        yule_walker(series, 1, phi);
        REQUIRE(phi.size() == 1u);
        // Should recover ≈ 0.8 (within ±0.15 given small noise)
        CHECK(phi[0] == doctest::Approx(0.8).epsilon(0.15));
    }

} // TEST_SUITE
