# Spec: Signal Kernels Forward Delivery Contract

## Objective

Bring Signal Kernels to the shared Project Telos public/developer delivery floor
while preserving its header-only C++23 algorithm behavior.

## Requirements

- [x] Keep README, USAGE, AGENTS, CMake build notes, tests, examples, and CI
  aligned.
- [x] Add a current changelog and implementation receipt.
- [x] Use native CMake on GitHub-hosted Windows runners with current checkout.
- [x] Normalize forward-facing punctuation so the public-surface scanner reports
  a clean boundary.

## Technical Approach

Use a documentation and CI-only patch. Existing C++ tests remain the behavioral
authority.

## Success Criteria

- [x] `cmake --build build --config Debug` passes.
- [x] `ctest --test-dir build -C Debug --output-on-failure` passes.
- [x] `python -m public_surface_sweeper . --workspace --json` reports `MATCH`.
- [x] `git diff --check` exits 0.

## Status: IMPLEMENTED
