# Changelog

## 1.1.0 (unreleased)

- From v1.1.0, code is licensed FSL-1.1-MIT. Earlier releases remain under MIT.
- `LICENSE` is the FSL-1.1-MIT text from fsl.software, with licensor Zain Dana Harper and copyright 2026.
- The project version in `CMakeLists.txt` moves from 1.0.0 to 1.1.0. No version was tagged before this change, so MIT covers every commit before it.
- `tests/third_party/doctest/` keeps doctest's own MIT licence, copyright Viktor Kirilov. Its `LICENSE.txt`, which the header names, was missing and is now added from doctest v2.4.11.
- No code changed.

## 2026-06-29 - Forward Delivery Contract

- Added `project-docs/specs/SPEC-signal-kernels-forward-delivery.md` as the
  implementation receipt for the delivery pass.
- Simplified CI to use native CMake on `windows-latest` plus current
  `actions/checkout`.
- Normalized forward-facing punctuation for public-surface scanner
  compatibility.
- Kept header algorithms, test coverage, CMake targets, and public API behavior
  unchanged.

## Current Status

- Runtime: header-only C++23 library.
- Surfaces: public headers, CMake INTERFACE target, tests, usage guide, and
  demo pipeline.
- Verification: CMake build plus CTest.
