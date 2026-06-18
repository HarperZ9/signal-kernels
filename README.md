# Signal Kernels

> Header-only C++23 signal/information-theory library: entropy, MI/transfer entropy, divergences, Granger, PELT, FFT, and forecasting.

[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
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

- `algorithms/entropy.hpp` — Shannon, Rényi, Tsallis, min, block, spectral, and
  permutation entropy.
- `algorithms/information.hpp` — mutual information, transfer entropy, and
  KL / Jensen-Shannon / Hellinger / Wasserstein divergences.
- `algorithms/causal.hpp` — Granger causality test.
- `algorithms/changepoint.hpp` — PELT (Pruned Exact Linear Time) change-point
  detection.
- `algorithms/forecast.hpp` — SARIMA and VAR time-series forecasting.
- `algorithms/curvature.hpp` — Forman-Ricci and Ollivier-Ricci graph curvature.
- `algorithms/_fft.hpp` — radix-2 Cooley-Tukey FFT (internal).
- `algorithms/_numeric.hpp` — Welford variance, log-sum-exp, small-matrix linear
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

## Building and testing

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## License

MIT. See [LICENSE](LICENSE).

---
**Zain Dana Harper** — small tools with explicit edges.
[Portfolio](https://harperz9.github.io) · [HarperZ9](https://github.com/HarperZ9)
<sub>Built with Claude Code; reviewed, tested, and owned by me.</sub>
