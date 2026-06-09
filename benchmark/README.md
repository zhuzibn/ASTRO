# ASTRO Reproducible Baseline Benchmark

## Purpose

This benchmark freezes the current ASTRO simulation behavior before
refactoring. It follows the audited `main.m` parameter and solver path while
replacing the non-reproducible `rng('shuffle')` atom distribution with a saved
fixed matrix.

This is a regression reference, not a proof or full validation of the
underlying physics.

## Run

Run commands from the repository root.

Quick `3 x 5` baseline (default):

```bash
matlab -batch "run('benchmark/run_baseline_benchmark.m')"
```

Current-size `20 x 30` baseline:

```bash
matlab -batch "benchmark_mode='current'; run('benchmark/run_baseline_benchmark.m')"
```

Both modes retain the audited physics and time settings. Quick mode changes
only `natomW` and `natomL`.

## Outputs

- `benchmark/baseline/input/baseline_atomtype.mat`: fixed current-size
  distribution.
- `benchmark/baseline/input/baseline_atomtype_quick.mat`: fixed quick
  distribution.
- `benchmark/baseline/output/baseline_result.mat`: explicit result schema for
  the most recently run mode.
- `benchmark/baseline/baseline_metadata.md`: environment, parameters, run
  status, warnings, and deviations.
- `benchmark/baseline/figures/atomtype_matrix.png`
- `benchmark/baseline/figures/initial_mz.png`
- `benchmark/baseline/figures/final_mz.png`
- `benchmark/baseline/figures/average_spin_vs_time.png`
- `benchmark/baseline/figures/spin_norm_histogram.png`

The runner does not depend on or overwrite root-level `final.mat`.

## Compare

From MATLAB:

```matlab
addpath('benchmark')
passed = compare_to_baseline('path/to/future_result.mat');
```

The comparison checks array sizes, exact atom distribution, final average
spin, full saved spin arrays, and final spin-norm statistics. Default
tolerances are:

```matlab
tol_avg_spin = 1e-8;
tol_spin_array = 1e-6;
tol_spin_norm = 1e-8;
```

These tolerances are conservative starting points for GPU runs. Record
observed platform differences before changing them.

## Fixed Distribution

`atomtype_ == 1` denotes RE/Gd and `atomtype_ == 0` denotes TM/FeCo. The
primary baseline uses `compositionn = 0.1`, matching the audited current code.
The presentation describes Co0.7Gd0.3, but the current code does not implement
that composition or distinguish Fe from Co sites. The benchmark deliberately
preserves current code behavior instead of changing it to match the
presentation.

The MAT-files store the generated matrix, dimensions, composition, seed, and
description. Future runs load the saved matrix; the seed is provenance rather
than the reproducibility mechanism.

## Plot Orientation

Spin matrices are indexed `(W, L, time)`. Figures map physical length `L` to
the x-axis and physical width `W` to the y-axis. The displayed matrix is
explicitly reversed with `flipud` so increasing physical W follows the
established plotting convention rather than MATLAB row display order.

The benchmark figures are self-contained and do not use
`plot_dat_proc/threedplott.m` or its missing `arrow3` dependency.

## Known Limitations

- GPU numerical results are not guaranteed to reproduce bit-for-bit across
  MATLAB, CUDA, driver, or GPU versions.
- Current ASTRO output has no explicit energy observable.
- The benchmark preserves current `compositionn = 0.1`; the presentation's
  Co0.7Gd0.3 description differs from the audited code.
- The existing saved-state/time-vector alignment may need future review.
- Width must be reversed when displaying matrix data in the physical
  orientation.
- The current solver path requires a compatible NVIDIA GPU and MATLAB GPU
  support. The benchmark does not change physics to provide a CPU fallback.
- Separate Gd and TM averages are included only as diagnostic quantities, not
  validated physics observables.
