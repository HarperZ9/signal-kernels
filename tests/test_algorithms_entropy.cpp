// =============================================================================
// Tests: algorithms — entropy module
//
// Known-value checks:
//   - Shannon uniform-8 = 3.0 bits
//   - Shannon 50/50 = 1.0 bit
//   - Rényi → Shannon as alpha → 1
//   - Tsallis → 0 for deterministic
//   - min_entropy = -log2(max_p)
//   - block_entropy uniform pairs = 1.0 bit
//   - spectral_entropy > 0 for non-trivial signal
//   - permutation_entropy > 0 for non-constant series
// =============================================================================

#include <doctest/doctest.h>
#include "algorithms/entropy.hpp"

#include <cmath>
#include <vector>

using namespace warden::algorithms;

TEST_SUITE("algorithms::entropy") {

    TEST_CASE("shannon uniform 8 symbols = 3.0 bits") {
        std::vector<double> p = {0.125, 0.125, 0.125, 0.125,
                                  0.125, 0.125, 0.125, 0.125};
        CHECK(shannon(p) == doctest::Approx(3.0).epsilon(1e-9));
    }

    TEST_CASE("shannon binary 50/50 = 1.0 bit") {
        std::vector<double> p = {0.5, 0.5};
        CHECK(shannon(p) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("shannon deterministic = 0") {
        std::vector<double> p = {1.0, 0.0, 0.0};
        CHECK(shannon(p) == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("shannon_from_counts matches shannon") {
        std::vector<uint64_t> c = {1, 1, 1, 1, 1, 1, 1, 1};
        CHECK(shannon_from_counts(c) == doctest::Approx(3.0).epsilon(1e-9));
    }

    TEST_CASE("shannon_from_bytes saturated (all 256 symbols) ≈ 8 bits") {
        std::vector<uint8_t> bytes(256);
        for (int i = 0; i < 256; ++i) bytes[static_cast<size_t>(i)] = static_cast<uint8_t>(i);
        // Exactly one of each byte => H = log2(256) = 8 bits
        CHECK(shannon_from_bytes(bytes) == doctest::Approx(8.0).epsilon(1e-9));
    }

    TEST_CASE("shannon_from_bytes all-same = 0") {
        std::vector<uint8_t> bytes(64, 0xAA);
        CHECK(shannon_from_bytes(bytes) == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("renyi alpha=1 matches shannon") {
        std::vector<double> p = {0.25, 0.25, 0.25, 0.25};
        double h_s = shannon(p);
        double h_r = renyi(p, 1.0);
        CHECK(h_r == doctest::Approx(h_s).epsilon(1e-9));
    }

    TEST_CASE("renyi alpha=2 collision entropy <= shannon") {
        std::vector<double> p = {0.5, 0.3, 0.2};
        double h2 = renyi(p, 2.0);
        double hs = shannon(p);
        CHECK(h2 <= hs + 1e-9);
    }

    TEST_CASE("tsallis q=1 matches shannon (epsilon check)") {
        std::vector<double> p = {0.25, 0.25, 0.25, 0.25};
        // q→1 limit, we pass q=1.0 which uses the shannon branch
        double t1 = tsallis(p, 1.0);
        double hs = shannon(p);
        CHECK(t1 == doctest::Approx(hs).epsilon(1e-6));
    }

    TEST_CASE("tsallis q=2 deterministic = 0") {
        std::vector<double> p = {1.0};
        CHECK(tsallis(p, 2.0) == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("min_entropy = -log2(max_p)") {
        std::vector<double> p = {0.5, 0.3, 0.2};
        CHECK(min_entropy(p) == doctest::Approx(-std::log2(0.5)).epsilon(1e-9));
    }

    TEST_CASE("min_entropy uniform n=4 = 2 bits") {
        std::vector<double> p = {0.25, 0.25, 0.25, 0.25};
        CHECK(min_entropy(p) == doctest::Approx(2.0).epsilon(1e-9));
    }

    TEST_CASE("block_entropy 2-byte pairs uniform = 1.0 bit") {
        // Two distinct blocks: {0x00,0x01} and {0x02,0x03}, 4 of each => H=1
        std::vector<uint8_t> seq = {0, 1, 0, 1, 0, 1, 0, 1,
                                     2, 3, 2, 3, 2, 3, 2, 3};
        double h = block_entropy(seq, 2);
        CHECK(h == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("block_entropy single block = 0") {
        std::vector<uint8_t> seq = {0xAB, 0xCD, 0xAB, 0xCD, 0xAB, 0xCD};
        double h = block_entropy(seq, 2);
        // Only one distinct 2-byte block
        CHECK(h == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("spectral_entropy non-zero for sine wave") {
        constexpr size_t N = 128;
        std::vector<double> sig(N);
        for (size_t i = 0; i < N; ++i)
            sig[i] = std::sin(2.0 * 3.14159265358979323846 * 8.0 * static_cast<double>(i) / N);
        double se = spectral_entropy(sig, 1.0);
        CHECK(se > 0.0);
    }

    TEST_CASE("spectral_entropy pure tone < white noise") {
        constexpr size_t N = 256;
        std::vector<double> tone(N), noise(N);
        for (size_t i = 0; i < N; ++i) {
            tone[i]  = std::sin(2.0 * 3.14159265358979323846 * 4.0 * static_cast<double>(i) / N);
            noise[i] = static_cast<double>((i * 1664525u + 1013904223u) & 0xFFFFFFu) / 0x1000000u - 0.5;
        }
        CHECK(spectral_entropy(tone, 1.0) < spectral_entropy(noise, 1.0));
    }

    TEST_CASE("permutation_entropy > 0 for non-constant series") {
        std::vector<double> s = {0.1, 0.5, 0.3, 0.9, 0.2, 0.7, 0.4, 0.8, 0.6};
        double pe = permutation_entropy(s, 3, 1);
        CHECK(pe > 0.0);
    }

    TEST_CASE("permutation_entropy monotone = 0") {
        // Strictly increasing => only one ordinal pattern => H = 0
        std::vector<double> s;
        for (int i = 0; i < 20; ++i) s.push_back(static_cast<double>(i));
        double pe = permutation_entropy(s, 3, 1);
        CHECK(pe == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("permutation_entropy bounded by log2(order!)") {
        std::vector<double> s = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3};
        double pe = permutation_entropy(s, 3, 1);
        double upper = std::log2(6.0); // 3! = 6
        CHECK(pe <= upper + 1e-9);
    }

} // TEST_SUITE
