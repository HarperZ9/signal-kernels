<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/HarperZ9/signal-kernels/main/docs/art/hero-dark.svg">
  <img src="https://raw.githubusercontent.com/HarperZ9/signal-kernels/main/docs/art/hero-light.svg" alt="signal-kernels: C++23 header-only entropy, causality, and forecasting kernels. 6 wavering traces run from the left and narrow into a bright core over a row of tick marks." width="100%">
</picture>

# signal-kernels

C++23 header-only entropy, causality, and forecasting kernels.

```
cmake -S . -B build -DSIGNAL_KERNELS_BUILD_TESTS=ON
```

[![version: 1.1.0](https://img.shields.io/badge/version-1.1.0-e6e1d6?style=flat-square&labelColor=1a1712)](https://github.com/HarperZ9/signal-kernels/releases/latest)
[![CI](https://github.com/HarperZ9/signal-kernels/actions/workflows/ci.yml/badge.svg)](https://github.com/HarperZ9/signal-kernels/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-FSL--1.1--MIT-e6e1d6?style=flat-square&labelColor=1a1712)](https://github.com/HarperZ9/signal-kernels/blob/main/LICENSE)
![C++23](https://img.shields.io/badge/language-C%2B%2B23-e6e1d6?style=flat-square&labelColor=1a1712)

Signal Kernels is a header-only C++23 library for scientific signal processing
and telemetry analytics. It includes entropy measures, divergence metrics,
causal tests, change-point detection, FFT helpers, graph curvature, and
forecasting primitives.

## See it work, step by step

The [animated explainer](https://harperz9.github.io/repo-explainers/signal-kernels.html)
walks through each module of the library on the inputs of examples/demo_pipeline.cpp: entropy, divergences, Granger causality, PELT change points, SARIMA and VAR forecasts, and graph curvature. Every value on it is output from this repository. Its
source is [docs/explainer/index.html](docs/explainer/index.html).

## Watch

No concept film fits this tool closely yet. The walkthrough below covers it in text, with real commands and output.

Video walkthrough: coming with the next release.

## Walkthrough

Install it, run it once, then use the main feature. Each command below is real, and so is its output.

1. **Get it and build.** Clone and build with CMake and MSVC; the bundled CMake file targets Windows x64.

   ```text
   $ git clone https://github.com/HarperZ9/signal-kernels && cd signal-kernels
   $ cmake -S . -B build -DSIGNAL_KERNELS_BUILD_TESTS=ON
   $ cmake --build build --config Debug
   ```

2. **Run the tests.** The test binary covers every header.

   ```text
   $ ctest --test-dir build -C Debug --output-on-failure
   ```

3. **Entropy.** The demo program `examples/demo_pipeline.cpp`, built in Release, printed these entropy values.

   ```text

   shannon(uniform-8)        = 3.000000 bits
   renyi(uniform-8, a=2)     = 3.000000 bits
   min_entropy(uniform-8)    = 3.000000 bits
   shannon_from_bytes(0..255) = 8.000000 bits
   permutation_entropy(o=3)  = 1.842371 bits
   ```

4. **Change points.** On a step series, PELT finds one change point at index 25.

   ```text
   input: step = 25 x 0.0, then 25 x 10.0
   pelt(L2) detected 1 change point(s):
     index=25  segment_cost=0.000000
   ```

## Why it matters

Large AI and research systems need reliable measurement kernels before a model
interprets noisy data. This repo provides public, testable primitives that can
feed receipt-backed scientific and operational workflows.

## Try it

```bash
cmake -S . -B build -DSIGNAL_KERNELS_BUILD_TESTS=ON
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

## What to test first

- Build the CMake project with tests enabled.
- Run CTest.
- Include `algorithms/entropy.hpp` or another module in a small C++ target.

## Current status

Header-only C++23 library with unit tests. The bundled CMake setup targets
Windows x64/MSVC, while the headers are standard-library-only.

## Existing technical notes

> Header-only C++23 signal/information-theory library: entropy, MI/transfer entropy, divergences, Granger, PELT, FFT, and forecasting.

[![license: FSL-1.1-MIT](https://img.shields.io/badge/license-FSL--1.1--MIT-blue.svg)](LICENSE)
![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)
[![CI](https://github.com/HarperZ9/signal-kernels/actions/workflows/ci.yml/badge.svg)](https://github.com/HarperZ9/signal-kernels/actions/workflows/ci.yml)
![header-only](https://img.shields.io/badge/header--only-C%2B%2B23-success.svg)
[![part of: AI-accountability toolkit](https://img.shields.io/badge/part_of-AI--accountability_toolkit-7a5cff.svg)](https://harperz9.github.io)

## Overview

`signal-kernels` implements a focused set of well-known, published analytics
primitives for entropy, forecasting, causal analysis, change-point detection,
graph curvature, and spectral analysis. It depends only on the C++ standard
library and is intended as a reusable foundation for telemetry analytics and
scientific signal processing.

## Modules

- `algorithms/entropy.hpp` -- Shannon, Rényi, Tsallis, min, block, spectral, and
  permutation entropy.
- `algorithms/information.hpp` -- mutual information, transfer entropy, and
  KL / Jensen-Shannon / Hellinger / Wasserstein divergences.
- `algorithms/causal.hpp` -- Granger causality test.
- `algorithms/changepoint.hpp` -- PELT (Pruned Exact Linear Time) change-point
  detection.
- `algorithms/forecast.hpp` -- SARIMA and VAR time-series forecasting.
- `algorithms/curvature.hpp` -- Forman-Ricci and Ollivier-Ricci graph curvature.
- `algorithms/_fft.hpp` -- radix-2 Cooley-Tukey FFT (internal).
- `algorithms/_numeric.hpp` -- Welford variance, log-sum-exp, small-matrix linear
  algebra, autocorrelation, Yule-Walker (internal).

## Why this is publishable

- No offensive-action primitives (no exploitation, command execution,
  persistence, lateral movement, or credential tooling).
- No secrets, credentials, keys, or operator-specific identifiers.
- Standalone build system and public API header set.
- Unit tests are included for all headers.

## Usage

Add `signal-kernels` to your CMake project and link the `signal-kernels`
INTERFACE target:

```cmake
add_subdirectory(signal-kernels)
target_link_libraries(your_target PRIVATE signal-kernels)
```

Then include the headers you need:

```cpp
#include "algorithms/entropy.hpp"

double h = signal_kernels::algorithms::shannon(probs);
```

See [USAGE.md](USAGE.md) for a full walkthrough -- the public function/class
list per header, worked examples with expected output, and build notes. A
single end-to-end program lives at
[`examples/demo_pipeline.cpp`](examples/demo_pipeline.cpp).

## Platform

The bundled `CMakeLists.txt` targets **Windows x64 / MSVC only** and stops with
a fatal error on other platforms. The library is header-only and uses only the
C++23 standard library, so the headers can be compiled with other conforming
toolchains if you bypass the bundled CMake configuration.

## Building and testing

```bash
cmake -S . -B build -DSIGNAL_KERNELS_BUILD_TESTS=ON
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Tests are gated on `SIGNAL_KERNELS_BUILD_TESTS` (defaults to `ON` only when
`signal-kernels` is the top-level project). They use
[doctest](https://github.com/doctest/doctest): a vendored copy at
`tests/third_party/doctest/doctest.h` is used if present, otherwise doctest
`v2.4.11` is fetched via `FetchContent`. This test-only dependency does not
affect consumers of the header-only library.

## License

From v1.1.0, code is licensed FSL-1.1-MIT. Earlier releases remain under MIT. FSL-1.1-MIT is the Functional Source License, Version 1.1, with MIT as the future licence: each release becomes available under MIT two years after it is made available. See [LICENSE](LICENSE). Every commit in this repository is by the author.

The vendored test framework `tests/third_party/doctest/doctest.h` is doctest 2.4.11 by Viktor Kirilov, under its own MIT licence ([tests/third_party/doctest/LICENSE.txt](tests/third_party/doctest/LICENSE.txt)). It is test-only and not relicensed.

---
**Zain Dana Harper** -- small tools with explicit edges.
[Portfolio](https://harperz9.github.io) · [HarperZ9](https://github.com/HarperZ9)
<sub>Built with Claude Code; reviewed, tested, and owned by me.</sub>

## For developers

Keep the public README, build notes, and examples aligned with current behavior. Before opening a PR or pushing a release, run the local native verification path.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

See [AGENTS.md](AGENTS.md) for the repo-specific operating boundary and
[CHANGELOG.md](CHANGELOG.md) for current delivery status.

---

Built by **[Zain Dana Harper](https://harperz9.github.io)** in Seattle: evidence-first tools that leave a re-checkable artifact behind. The full workbench is at [Project Telos](https://harperz9.github.io).
