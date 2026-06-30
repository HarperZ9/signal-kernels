// Best-effort demo -- not runtime-verified by author.
// =============================================================================
// examples/demo_pipeline.cpp
//
// End-to-end tour of the signal-kernels public API. Every call below uses an
// existing public function or type from include/algorithms/. Output is printed
// to stdout; nothing is read from disk or the network.
//
// Build (header-only -- just put include/ on the search path), e.g. with MSVC:
//   cl /std:c++latest /EHsc /I ..\include demo_pipeline.cpp
//
// Or add this file as an executable in your own CMake project that links the
// `signal-kernels` INTERFACE target.
// =============================================================================

#include "algorithms/entropy.hpp"
#include "algorithms/information.hpp"
#include "algorithms/causal.hpp"
#include "algorithms/changepoint.hpp"
#include "algorithms/forecast.hpp"
#include "algorithms/curvature.hpp"

#include <cstdint>
#include <cstdio>
#include <vector>

using namespace signal_kernels::algorithms;

// Deterministic LCG so the demo is reproducible without <random> seeding noise.
static double next_unit(uint64_t& state) {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>((state >> 11) & 0xFFFFFF) / 16777216.0; // [0,1)
}

int main() {
    // -------------------------------------------------------------------------
    // 1. Entropy
    // -------------------------------------------------------------------------
    std::printf("== entropy ==\n");
    std::vector<double> uniform8(8, 0.125);
    std::printf("shannon(uniform-8)       = %.6f bits (expect 3.0)\n",
                shannon(uniform8));
    std::printf("renyi(uniform-8, a=2)    = %.6f bits\n", renyi(uniform8, 2.0));
    std::printf("min_entropy(uniform-8)   = %.6f bits\n", min_entropy(uniform8));

    std::vector<uint8_t> bytes(256);
    for (int i = 0; i < 256; ++i) bytes[static_cast<size_t>(i)] =
        static_cast<uint8_t>(i);
    std::printf("shannon_from_bytes(0..255)= %.6f bits (expect 8.0)\n",
                shannon_from_bytes(bytes));

    std::vector<double> wiggly = {0.1, 0.5, 0.3, 0.9, 0.2, 0.7, 0.4, 0.8, 0.6};
    std::printf("permutation_entropy(o=3) = %.6f bits\n",
                permutation_entropy(wiggly, 3, 1));

    // -------------------------------------------------------------------------
    // 2. Information-theoretic divergences
    // -------------------------------------------------------------------------
    std::printf("\n== information ==\n");
    std::vector<double> p = {0.5, 0.5};
    std::vector<double> q = {0.9, 0.1};
    std::printf("kl_divergence(p||q)      = %.6f bits\n", kl_divergence(p, q));
    std::printf("js_divergence(p,q)       = %.6f bits\n", js_divergence(p, q));
    std::printf("hellinger(p,q)           = %.6f\n",      hellinger(p, q));

    std::vector<double> u = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> v = {2.0, 3.0, 4.0, 5.0};
    std::printf("wasserstein_1d(u,v)      = %.6f\n", wasserstein_1d(u, v));

    std::vector<double> xa = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3};
    std::printf("mutual_information(x,x)  = %.6f bits\n",
                mutual_information(xa, xa, 4));

    // -------------------------------------------------------------------------
    // 3. Granger causality:  build y_t = 0.8 * x_{t-1} + noise
    // -------------------------------------------------------------------------
    std::printf("\n== causal ==\n");
    constexpr int N = 200;
    uint64_t rng = 0xDEADBEEF12345678ull;
    std::vector<double> gx(N), gy(N, 0.0);
    for (int i = 0; i < N; ++i) gx[static_cast<size_t>(i)] = next_unit(rng);
    for (int i = 1; i < N; ++i)
        gy[static_cast<size_t>(i)] =
            0.8 * gx[static_cast<size_t>(i - 1)] + 0.01 * next_unit(rng);

    GrangerResult gr = granger_causality(gx, gy, 3);
    std::printf("granger x->y: f_stat=%.4f  p_value=%.4f  optimal_lag=%d\n",
                gr.f_stat, gr.p_value, gr.optimal_lag);

    // -------------------------------------------------------------------------
    // 4. PELT change-point detection on a step signal (step at index 25)
    // -------------------------------------------------------------------------
    std::printf("\n== changepoint ==\n");
    std::vector<double> step(50, 0.0);
    for (size_t i = 25; i < 50; ++i) step[i] = 10.0;

    auto cps = pelt(step, cost_l2, bic_penalty(step.size()), 2);
    std::printf("pelt(L2) detected %zu change point(s):\n", cps.size());
    for (const auto& cp : cps)
        std::printf("  index=%zu  segment_cost=%.6f\n", cp.index, cp.cost);

    // -------------------------------------------------------------------------
    // 5. Forecasting: SARIMA(1,0,0) and a 2-series VAR(1)
    // -------------------------------------------------------------------------
    std::printf("\n== forecast ==\n");
    std::vector<double> series(200, 0.0);
    for (size_t i = 1; i < series.size(); ++i)
        series[i] = 0.7 * series[i - 1] + 0.01 * static_cast<double>(i % 5);

    SARIMA sarima(1, 0, 0, 0, 0, 0, 1);
    sarima.fit(series);
    std::printf("SARIMA fitted=%d  ar[0]=%.4f  sigma2=%.6f\n",
                sarima.is_fitted() ? 1 : 0,
                sarima.ar_params().empty() ? 0.0 : sarima.ar_params()[0],
                sarima.sigma2());
    auto sf = sarima.forecast(5);
    std::printf("SARIMA forecast(5) -> %zu values:", sf.size());
    for (double f : sf) std::printf(" %.4f", f);
    std::printf("\n");

    std::vector<double> a(100), b(100);
    uint64_t rng2 = 0x1234567890ABCDEFull;
    for (size_t i = 1; i < 100; ++i) {
        a[i] = 0.6 * a[i - 1] + 0.02 * next_unit(rng2);
        b[i] = 0.4 * b[i - 1] + 0.02 * next_unit(rng2);
    }
    VAR var(1);
    var.fit({a, b});
    auto vf = var.forecast(3);
    std::printf("VAR(1) fitted=%d  forecast shape = [%zu][%zu]\n",
                var.is_fitted() ? 1 : 0,
                vf.size(), vf.empty() ? 0 : vf[0].size());

    // -------------------------------------------------------------------------
    // 6. Graph curvature on a 4-cycle C4
    // -------------------------------------------------------------------------
    std::printf("\n== curvature ==\n");
    Graph g;
    g.add_edge(0, 1, 1.0);
    g.add_edge(1, 2, 1.0);
    g.add_edge(2, 3, 1.0);
    g.add_edge(3, 0, 1.0);
    std::printf("forman_ricci(0,1)        = %.6f\n", forman_ricci(g, 0, 1));
    std::printf("ollivier_ricci(0,1,a=0.5)= %.6f\n", ollivier_ricci(g, 0, 1, 0.5));
    std::printf("shortest_path(0->2)      = %.6f\n", g.shortest_path(0, 2));

    std::printf("\nDone.\n");
    return 0;
}
