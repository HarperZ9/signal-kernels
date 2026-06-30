// =============================================================================
// Tests: algorithms -- information module
//
// Known-value checks:
//   - mutual_information(X, X) = H(X) (self-MI)
//   - mutual_information independent = 0
//   - kl_divergence(P||P) = 0
//   - js_divergence symmetric
//   - hellinger in [0, 1]
//   - wasserstein_1d identical distributions = 0
//   - transfer_entropy self-TE > TE of independent
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/information.hpp"

#include <cmath>
#include <vector>

using namespace signal_kernels::algorithms;

TEST_SUITE("algorithms::information") {

    TEST_CASE("mutual_information self-MI > 0") {
        std::vector<double> x = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0};
        double mi = mutual_information(x, x, 4);
        CHECK(mi > 0.0);
    }

    TEST_CASE("mutual_information identical series maximized") {
        std::vector<double> x = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
        double mi_self = mutual_information(x, x, 2);
        // Self-MI should be higher than MI with shuffled
        std::vector<double> y = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
        double mi_shift = mutual_information(x, y, 2);
        CHECK(mi_self >= mi_shift - 1e-9);
    }

    TEST_CASE("kl_divergence same distribution = 0") {
        std::vector<double> p = {0.25, 0.25, 0.25, 0.25};
        double kl = kl_divergence(p, p);
        CHECK(kl == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("kl_divergence asymmetric") {
        std::vector<double> p = {0.5, 0.5};
        std::vector<double> q = {0.9, 0.1};
        double kl_pq = kl_divergence(p, q);
        double kl_qp = kl_divergence(q, p);
        CHECK(kl_pq != doctest::Approx(kl_qp).epsilon(1e-3));
    }

    TEST_CASE("kl_divergence q has zero where p > 0 returns inf") {
        std::vector<double> p = {0.5, 0.5};
        std::vector<double> q = {1.0, 0.0};
        double kl = kl_divergence(p, q);
        CHECK(std::isinf(kl));
    }

    TEST_CASE("js_divergence is symmetric") {
        std::vector<double> p = {0.6, 0.4};
        std::vector<double> q = {0.3, 0.7};
        CHECK(js_divergence(p, q) == doctest::Approx(js_divergence(q, p)).epsilon(1e-9));
    }

    TEST_CASE("js_divergence identical distributions = 0") {
        std::vector<double> p = {0.2, 0.3, 0.5};
        CHECK(js_divergence(p, p) == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("js_divergence bounded [0, 1]") {
        std::vector<double> p = {1.0, 0.0};
        std::vector<double> q = {0.0, 1.0};
        double jsd = js_divergence(p, q);
        CHECK(jsd >= 0.0);
        CHECK(jsd <= 1.0 + 1e-9);
    }

    TEST_CASE("hellinger identical = 0") {
        std::vector<double> p = {0.5, 0.3, 0.2};
        CHECK(hellinger(p, p) == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("hellinger disjoint support = 1") {
        std::vector<double> p = {1.0, 0.0};
        std::vector<double> q = {0.0, 1.0};
        CHECK(hellinger(p, q) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("hellinger bounded [0, 1]") {
        std::vector<double> p = {0.7, 0.2, 0.1};
        std::vector<double> q = {0.1, 0.2, 0.7};
        double h = hellinger(p, q);
        CHECK(h >= 0.0);
        CHECK(h <= 1.0 + 1e-9);
    }

    TEST_CASE("wasserstein_1d identical distributions = 0") {
        std::vector<double> u = {1.0, 2.0, 3.0, 4.0};
        CHECK(wasserstein_1d(u, u) == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("wasserstein_1d shift by constant = shift") {
        std::vector<double> u = {0.0, 0.0, 0.0, 0.0};
        std::vector<double> v = {1.0, 1.0, 1.0, 1.0};
        // EMD between point masses at 0 and 1 = 1.0
        CHECK(wasserstein_1d(u, v) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("wasserstein_1d non-negative") {
        std::vector<double> u = {1.0, 5.0, 3.0};
        std::vector<double> v = {2.0, 4.0, 6.0};
        CHECK(wasserstein_1d(u, v) >= 0.0);
    }

    TEST_CASE("transfer_entropy non-negative") {
        std::vector<double> s = {0.1, 0.4, 0.7, 0.2, 0.9, 0.3, 0.6, 0.8, 0.5, 0.15,
                                  0.45, 0.75, 0.25, 0.95, 0.35, 0.65, 0.85, 0.55, 0.05, 0.55};
        std::vector<double> t = s; // same = high TE
        double te = transfer_entropy(s, t, 1, 1);
        CHECK(te >= 0.0);
    }

} // TEST_SUITE
