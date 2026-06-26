# Changelog

## Source Code Changes

### 2026-06-25 - Extract deterministic spin initialization

- Added `astro_validate_atom_distribution.m` to validate deterministic RE/TM
  atom matrices for shape, allowed values, and optional composition.
- Added `astro_initial_spin_state.m` to construct deterministic initial spin
  arrays from explicit inputs while preserving the existing RE/Gd, TM/FeCo,
  and domain-wall sign conventions.
- Updated `systemgeneration.m` and `benchmark/run_baseline_benchmark.m` to use
  the shared initializer while keeping random, fixed-input, and restart policy
  in their existing callers.
- Extended smoke checks with exact-output tests against the previous
  deterministic initialization formula.

### 2026-06-25 - Share production default configuration

- Added `astro_default_config.m` as the single source for ASTRO production
  physics, lattice, solver, torque, and time defaults.
- Updated `main.m` and `benchmark/run_baseline_benchmark.m` to consume the
  shared defaults, with benchmark overrides limited to lattice size, fixed
  atom distribution, and candidate output paths.
- Extended smoke checks to assert benchmark results still match the shared
  default configuration.

### 2026-06-25 - Reject unsupported solver modes

- Added an `ASTRO:UnsupportedSolver` guard in `main.m` and
  `integrate_llg.m` so only `rk4=1` reaches integration.
- Added smoke-test checks that `rk4=0` and `rk4=2` fail with the documented
  unsupported-solver message.
- Documented the RK4-only solver policy in `README.md`, `AGENTS.md`, and
  `benchmark/README.md`.

### 2026-06-25 - Add deterministic invariant smoke checks

- Added `tests/run_smoke_tests.m`, a MATLAB batch smoke entry point that runs
  the quick benchmark into a temporary candidate directory and exits nonzero
  if deterministic atom counts, trajectory dimensions, finite spin values,
  spin normalization, or saved-time invariants fail.
- Documented the smoke command in `README.md` and `AGENTS.md`.

### 2026-06-23 - Align saved spin states with saved times

- Updated the RK4 trajectory assembly in `integrate_llg.m` so the saved
  states include the initial state at `t=0`, include the final state at
  `runtime`, and skip duplicated carried states at chunk boundaries.
- Updated `main.m` and `benchmark/run_baseline_benchmark.m` to build the
  saved-time vector directly from `runtime`, `tstep`, and `savetstep`, making
  `numel(t)` match the saved spin trajectory depth.

### 2026-06-23 - Add non-destructive benchmark output mode

- Added optional `benchmark_output_root` handling to
  `benchmark/run_baseline_benchmark.m` so candidate benchmark result files,
  figures, and metadata can be written outside `benchmark/baseline/`.
- Documented the candidate output workflow in `benchmark/README.md`.

## Error Logs

### 2026-06-25 - Unsupported solver modes reached broken integration branches

- Symptom: selecting `rk4=0` or `rk4=2` could enter unvalidated Heun or
  predictor-corrector integration branches and fail later with misleading
  undefined-variable or obsolete-signature errors.
- Root cause: the production path only validated non-RK4 mode for
  spin-torque-driven runs, while unsupported solver branches remained
  reachable.
- Fix: added a shared `ASTRO:UnsupportedSolver` rejection before integration
  setup and covered `rk4=0` and `rk4=2` with deterministic smoke checks.
- Prevention: keep unsupported solver selections covered by negative tests
  until those solvers receive explicit validation and regression coverage.
