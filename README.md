# Warden Analytics Algorithms

## Public surface summary

`warden-algorithms` publishes a small, header-only subset of WARDEN analytics primitives for entropy, forecasting, causal and curvature analytics. It is intended as a reusable C++ foundation for telemetry analytics and scientific processing.

## Why this is publishable

- No offensive action primitives (no exploitation, command execution, persistence, lateral movement, or credential tooling).
- No secrets, credentials, keys, or operator-specific identifiers.
- Standalone build system and public API header set.
- Unit tests are included for all included headers.

## Modules exported

- `include/algorithms/*` and support headers.
- C++ tests in `tests/` that exercise the API contract.

## Public leaves and gates

- Test gate: `tests/test_algorithms_*.cpp` present.
- License gate: `LICENSE` present (MIT).
- Secret gate: secrets intentionally absent; keys, tokens, `.env`, credentials are not included.
- Claim gate: `PUBLIC-CLAIM.md` + `PUBLIC-DISCLAIMER.md` present.

## Usage

Add `warden-algorithms` to your CMake project and link to `we-algorithms` target.

```bash
cmake --build .
```

## Scope policy

Keep integrations defensive and analytical. Do not use this module for offensive signal abuse.
