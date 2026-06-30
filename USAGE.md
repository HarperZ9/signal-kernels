# Usage

`signal-kernels` is a **header-only C++23** library. There is no compiled
library to install and no command-line tool -- you consume it by adding its
`include/` directory to your build and `#include`-ing the headers you need.

All public functions and types live in the namespace
`signal_kernels::algorithms`.

> Expected outputs below marked **(illustrative)** were hand-derived from the
> algorithm definitions and the repository's unit tests; they were not produced
> by running the snippets through a compiler.

---

## Install / consume

### As a CMake subproject (recommended)

```cmake
add_subdirectory(signal-kernels)
target_link_libraries(your_target PRIVATE signal-kernels)
```

The `signal-kernels` target is an `INTERFACE` library: linking it just adds the
`include/` path and requests `cxx_std_23`. Then include headers directly:

```cpp
#include "algorithms/entropy.hpp"
```

### Without CMake

Add `include/` to your compiler's header search path and compile with C++23:

```bash
# MSVC (the supported toolchain -- see "Platform" below)
cl /std:c++latest /EHsc /I path\to\signal-kernels\include your_app.cpp
```

### Platform

The bundled `CMakeLists.txt` enforces **Windows x64 only** (`cmake` errors out
on other platforms) and adds MSVC flags (`/W4 /WX /permissive- /Zc:__cplusplus`).
The headers themselves are standard C++23 and depend only on the standard
library, so they can be compiled with other conforming toolchains if you bypass
the bundled CMake configuration.

---

## Headers and what they export

| Header | Public surface |
| --- | --- |
| `algorithms/entropy.hpp` | `shannon`, `shannon_from_counts`, `shannon_from_bytes`, `renyi`, `tsallis`, `min_entropy`, `block_entropy`, `spectral_entropy`, `permutation_entropy` |
| `algorithms/information.hpp` | `mutual_information`, `transfer_entropy`, `kl_divergence`, `js_divergence`, `hellinger`, `wasserstein_1d` |
| `algorithms/causal.hpp` | `granger_causality`, struct `GrangerResult { f_stat, p_value, optimal_lag }` |
| `algorithms/changepoint.hpp` | `pelt`, cost functors `cost_l2` / `cost_l1` / `cost_poisson` (and `cost_*_make(series)` factories), `bic_penalty`, `aic_penalty`, struct `Changepoint`, type alias `CostFn` |
| `algorithms/forecast.hpp` | classes `SARIMA`, `VAR` |
| `algorithms/curvature.hpp` | class `Graph`, `forman_ricci`, `ollivier_ricci`, struct `EdgeCurvature`, type alias `NodeId` (`uint32_t`) |

`algorithms/_fft.hpp` and `algorithms/_numeric.hpp` are internal helpers
(`signal_kernels::algorithms::detail` for the FFT). They are not part of the
stable public surface, though they compile and are unit-tested.

Most functions take `std::span<const double>` (or `std::span<const uint8_t>` for
the byte/sequence variants) so they bind to `std::vector`, `std::array`, or raw
buffers without copying. Entropy/divergence results are in **bits** (log base 2).

---

## Worked examples

### 1. Shannon entropy and a divergence

```cpp
#include "algorithms/entropy.hpp"
#include "algorithms/information.hpp"
#include <vector>
#include <print>

int main() {
    using namespace signal_kernels::algorithms;

    std::vector<double> uniform8(8, 0.125);     // 8 equally-likely symbols
    std::vector<double> p = {0.5, 0.5};
    std::vector<double> q = {0.9, 0.1};

    std::println("H(uniform8) = {} bits", shannon(uniform8)); // 3.0
    std::println("KL(p||q)    = {} bits", kl_divergence(p, q));
    std::println("JS(p,q)     = {} bits", js_divergence(p, q));
    return 0;
}
```

Expected output **(illustrative)**:

```
H(uniform8) = 3 bits
KL(p||q)    = 0.736966 bits
JS(p,q)     = 0.152699 bits
```

(`H` of 8 equiprobable symbols is exactly `log2(8) = 3`. `KL(p||q)` for
`p=[.5,.5]`, `q=[.9,.1]` is `0.5*log2(.5/.9) + 0.5*log2(.5/.1)`.)

---

### 2. Entropy directly from bytes

`shannon_from_bytes` builds the 256-symbol histogram for you, so you can measure
the byte-level entropy of any buffer (0 bits = constant, up to 8 bits = every
byte equally likely).

```cpp
#include "algorithms/entropy.hpp"
#include <vector>
#include <print>

int main() {
    using namespace signal_kernels::algorithms;

    std::vector<uint8_t> all_same(64, 0xAA);
    std::vector<uint8_t> all_distinct(256);
    for (int i = 0; i < 256; ++i) all_distinct[i] = static_cast<uint8_t>(i);

    std::println("constant buffer : {} bits", shannon_from_bytes(all_same));     // 0
    std::println("every byte once : {} bits", shannon_from_bytes(all_distinct)); // 8
    return 0;
}
```

Expected output **(illustrative)**:

```
constant buffer : 0 bits
every byte once : 8 bits
```

---

### 3. PELT change-point detection on a step signal

`pelt(series, cost_fn, penalty, min_size)` returns the detected change points
(each `Changepoint.index` is the first index of a new segment; index 0 and `n`
are not reported). Use a prebuilt cost functor (`cost_l2`, `cost_l1`,
`cost_poisson`) or the `cost_*_make(series)` factory for repeated calls on the
same series.

```cpp
#include "algorithms/changepoint.hpp"
#include <vector>
#include <print>

int main() {
    using namespace signal_kernels::algorithms;

    // 25 samples at 0.0, then 25 samples at 10.0 -- a step at index 25.
    std::vector<double> series(50, 0.0);
    for (size_t i = 25; i < 50; ++i) series[i] = 10.0;

    auto cps = pelt(series, cost_l2, bic_penalty(series.size()), 2);

    std::println("detected {} change point(s):", cps.size());
    for (const auto& cp : cps)
        std::println("  index={}  segment_cost={}", cp.index, cp.cost);
    return 0;
}
```

Expected output **(illustrative)**:

```
detected 1 change point(s):
  index=25  segment_cost=0
```

(The `test_algorithms_changepoint.cpp` "step function" case asserts the change
point lands within index 25 ±1.)

---

### 4. SARIMA forecast

`SARIMA(p, d, q, P, D, Q, s)` -- call `fit(series)` then `forecast(horizon)`.
With seasonal orders off, use `s = 1`. `forecast` always returns a vector of
length `horizon`; an unfitted model returns zeros.

```cpp
#include "algorithms/forecast.hpp"
#include <vector>
#include <print>

int main() {
    using namespace signal_kernels::algorithms;

    // Toy AR(1)-like series.
    std::vector<double> series(200, 0.0);
    for (size_t i = 1; i < series.size(); ++i)
        series[i] = 0.7 * series[i - 1] + 0.01 * static_cast<double>(i % 5);

    SARIMA model(/*p*/1, /*d*/0, /*q*/0, /*P*/0, /*D*/0, /*Q*/0, /*s*/1);
    model.fit(series);

    std::println("fitted    : {}", model.is_fitted());
    std::println("AR coef   : {}", model.ar_params().at(0));
    std::println("sigma^2   : {}", model.sigma2());

    auto fcst = model.forecast(5);
    std::println("forecast horizon = {}", fcst.size()); // 5
    return 0;
}
```

Expected output **(illustrative, values depend on the series)**:

```
fitted    : true
AR coef   : 0.69...
sigma^2   : 0.00...
forecast horizon = 5
```

(The forecast unit test only asserts `forecast(h).size() == h` and finiteness;
the recovered AR coefficient is checked to within ±10% of the generating phi.)

---

## Building and running the tests

```bash
cmake -S . -B build -DSIGNAL_KERNELS_BUILD_TESTS=ON
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Tests use [doctest](https://github.com/doctest/doctest). CMake uses a vendored
copy at `tests/third_party/doctest/doctest.h` if present, otherwise it fetches
doctest `v2.4.11` via `FetchContent` (so the test build -- unlike the library
itself -- requires network access on first configure unless doctest is vendored).
The test targets are gated on `SIGNAL_KERNELS_BUILD_TESTS`, which defaults to
`ON` only when `signal-kernels` is the top-level project.

See `examples/demo_pipeline.cpp` for a single end-to-end program that exercises
the entropy, information, change-point, forecasting, and curvature modules.
