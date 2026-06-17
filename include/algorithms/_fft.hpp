// Radix-2 Cooley-Tukey FFT
// =============================================================================
// algorithms/_fft.hpp -- Radix-2 Cooley-Tukey FFT (internal)
//
// Constraints:
//   - Input size must be a power of two, 1 ≤ N ≤ 65536 (2^16).
//   - Complex<double> I/O via std::complex<double>.
//   - Forward transform: exp(-2πi k n / N).
//   - Inverse transform: divide by N after forward on conjugate.
//
// Namespace: signal_kernels::algorithms::detail
// =============================================================================

#pragma once
#ifndef SIGNAL_KERNELS_FFT_HPP
#define SIGNAL_KERNELS_FFT_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace signal_kernels::algorithms::detail {

// Returns true iff n is a power of two (n > 0).
[[nodiscard]] constexpr bool is_power_of_two(size_t n) noexcept {
    return n > 0 && (n & (n - 1)) == 0;
}

// ---------------------------------------------------------------------------
// fft_inplace — in-place radix-2 DIT Cooley-Tukey FFT
//
// dir = -1  => forward  (standard DFT definition: exp(-2πi·k·n/N))
// dir = +1  => inverse  (caller must divide by N afterwards)
// ---------------------------------------------------------------------------

inline void fft_inplace(std::vector<std::complex<double>>& a, int dir) {
    const size_t N = a.size();
    if (N <= 1) return;
    if (!is_power_of_two(N))
        throw std::invalid_argument("fft_inplace: size must be a power of two");
    if (N > 65536)
        throw std::invalid_argument("fft_inplace: size exceeds 2^16 limit");

    // Bit-reversal permutation
    for (size_t i = 1, j = 0; i < N; ++i) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // Butterfly passes
    for (size_t len = 2; len <= N; len <<= 1) {
        double ang = dir * 2.0 * std::numbers::pi / static_cast<double>(len);
        std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < N; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (size_t j = 0; j < len / 2; ++j) {
                auto u = a[i + j];
                auto v = a[i + j + len / 2] * w;
                a[i + j]           = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// fft  — forward DFT (copy-in, returns complex spectrum)
// ifft — inverse DFT (returns time-domain signal, divided by N)
// magnitude_spectrum — returns |X[k]| for k = 0..N/2
// power_spectrum     — returns |X[k]|^2 / N (one-sided, DC + N/2 bins)
// ---------------------------------------------------------------------------

[[nodiscard]] inline std::vector<std::complex<double>>
fft(std::vector<std::complex<double>> x) {
    fft_inplace(x, -1);
    return x;
}

[[nodiscard]] inline std::vector<std::complex<double>>
ifft(std::vector<std::complex<double>> x) {
    fft_inplace(x, +1);
    const double inv = 1.0 / static_cast<double>(x.size());
    for (auto& v : x) v *= inv;
    return x;
}

[[nodiscard]] inline std::vector<double>
magnitude_spectrum(const std::vector<double>& real_signal) {
    const size_t n = real_signal.size();
    if (n == 0) return {};
    if (!is_power_of_two(n))
        throw std::invalid_argument("magnitude_spectrum: size must be power of two");

    std::vector<std::complex<double>> cx(n);
    for (size_t i = 0; i < n; ++i) cx[i] = {real_signal[i], 0.0};
    fft_inplace(cx, -1);

    // One-sided: DC + N/2 bins
    const size_t half = n / 2 + 1;
    std::vector<double> mag(half);
    for (size_t k = 0; k < half; ++k)
        mag[k] = std::abs(cx[k]);
    return mag;
}

[[nodiscard]] inline std::vector<double>
power_spectrum(const std::vector<double>& real_signal) {
    auto mag = magnitude_spectrum(real_signal);
    const double inv_n = 1.0 / static_cast<double>(real_signal.size());
    for (auto& m : mag) m = m * m * inv_n;
    return mag;
}

} // namespace signal_kernels::algorithms::detail

#endif // SIGNAL_KERNELS_FFT_HPP
