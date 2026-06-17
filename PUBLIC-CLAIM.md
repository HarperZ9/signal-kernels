# Public claim (signal-kernels)

- `signal-kernels` is a standalone, header-only C++ analytics library that
  implements well-known, published signal-processing and information-theory
  methods (entropy, mutual information, change-point detection, forecasting,
  graph curvature, and FFT-based spectral analysis).
- Scope is limited to statistical and mathematical utility components in
  `include/algorithms/`. There are no external dependencies beyond the C++
  standard library.
- The exported surface is a standalone C++23 CMake target with focused unit
  tests for every header.
- No credentials, operator notes, infrastructure targets, or private workflow
  material are present.
