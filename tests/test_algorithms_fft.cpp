// =============================================================================
// Tests: algorithms — _fft module
//
// Known-value checks:
//   - DC input => FFT output[0] = N * value, all others ≈ 0
//   - Single-tone => peak at correct bin
//   - Round-trip FFT then IFFT recovers signal
//   - Power spectrum peak at correct frequency bin
//   - Non-power-of-two throws
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/_fft.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

using namespace warden::algorithms::detail;

TEST_SUITE("algorithms::_fft") {

    TEST_CASE("DC signal: FFT[0] = N, rest ≈ 0") {
        constexpr size_t N = 8;
        std::vector<std::complex<double>> x(N, {1.0, 0.0});
        fft_inplace(x, -1);
        CHECK(std::abs(x[0]) == doctest::Approx(N).epsilon(1e-9));
        for (size_t k = 1; k < N; ++k)
            CHECK(std::abs(x[k]) < 1e-9);
    }

    TEST_CASE("single tone: peak at bin k") {
        constexpr size_t N = 64;
        constexpr size_t K = 7; // frequency bin
        std::vector<std::complex<double>> x(N);
        for (size_t n = 0; n < N; ++n) {
            double phase = 2.0 * std::numbers::pi * K * n / N;
            x[n] = {std::cos(phase), -std::sin(phase)};
        }
        fft_inplace(x, -1);
        double peak = std::abs(x[K]);
        CHECK(peak == doctest::Approx(N).epsilon(1e-6));
        for (size_t k = 0; k < N; ++k)
            if (k != K) CHECK(std::abs(x[k]) < 1e-6);
    }

    TEST_CASE("round-trip FFT then IFFT recovers signal") {
        constexpr size_t N = 16;
        std::vector<std::complex<double>> orig(N), x(N);
        for (size_t i = 0; i < N; ++i) orig[i] = x[i] = {std::sin(static_cast<double>(i)), 0.0};
        fft_inplace(x, -1);
        fft_inplace(x, +1);
        for (auto& v : x) v /= static_cast<double>(N);
        for (size_t i = 0; i < N; ++i)
            CHECK(x[i].real() == doctest::Approx(orig[i].real()).epsilon(1e-9));
    }

    TEST_CASE("magnitude_spectrum peak at correct bin") {
        constexpr size_t N = 128;
        constexpr size_t K = 5;
        std::vector<double> sig(N);
        for (size_t i = 0; i < N; ++i)
            sig[i] = std::cos(2.0 * std::numbers::pi * K * i / N);
        auto mag = magnitude_spectrum(sig);
        // Peak should be at bin K (one-sided)
        size_t peak_bin = 0;
        double peak_val = 0.0;
        for (size_t k = 0; k < mag.size(); ++k)
            if (mag[k] > peak_val) { peak_val = mag[k]; peak_bin = k; }
        CHECK(peak_bin == K);
    }

    TEST_CASE("power_spectrum sums are consistent") {
        constexpr size_t N = 32;
        std::vector<double> sig(N, 1.0);
        auto pwr = power_spectrum(sig);
        CHECK(!pwr.empty());
        CHECK(pwr[0] > 0.0); // DC power
    }

    TEST_CASE("non-power-of-two throws") {
        std::vector<std::complex<double>> x(7, {1.0, 0.0});
        CHECK_THROWS_AS(fft_inplace(x, -1), std::invalid_argument);
    }

    TEST_CASE("size > 2^16 throws") {
        std::vector<std::complex<double>> x(131072, {0.0, 0.0}); // 2^17
        CHECK_THROWS_AS(fft_inplace(x, -1), std::invalid_argument);
    }

    TEST_CASE("size 1 is a no-op") {
        std::vector<std::complex<double>> x = {{3.14, 0.0}};
        auto orig = x;
        fft_inplace(x, -1);
        CHECK(x[0].real() == doctest::Approx(orig[0].real()).epsilon(1e-12));
    }

} // TEST_SUITE
