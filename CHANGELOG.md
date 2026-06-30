# Changelog

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
