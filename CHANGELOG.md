# Changelog

## Source Code Changes

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

- None.
