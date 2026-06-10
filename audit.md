# ASTRO Repository Audit

Audit date: 2026-06-10

Scope: phase 1 repository audit. No MATLAB source, benchmark artifact, or
existing documentation was modified.

## Executive Summary

ASTRO is a MATLAB/GPU atomistic spin-dynamics model for a 2D ferrimagnetic
Gd-FeCo lattice. The active simulation is a script-driven pipeline:

```text
main.m
  -> constantfile.m
  -> systemgeneration.m
       -> distrib.m
  -> integrate_llg.m
       -> field_calc.m
            -> dipole_calc.m
       -> atomgpurk4.m / atomgpu.m / atomgpupc4.m
  -> final.mat
```

The default and benchmarked path uses RK4 on an NVIDIA GPU with exchange,
uniaxial anisotropy, DMI, and damping-like SOT enabled. Thermal and dipole
fields, STT, fixed edges, and external field are disabled or zero.

The repository has no build system, unit-test framework, CI configuration, or
automated lint command. It does have a reproducible regression benchmark with
fixed atom distributions, quick and current-size modes, explicit output
schema, figures, environment metadata, and a comparison function.

The checked-in current-size baseline was successfully generated on
2026-06-09. It stores a `20 x 30 x 100` trajectory, took 149.875465 seconds,
and has a final maximum spin-norm error of approximately `2.22e-16`.

## Repository State

- Branch: `master`
- HEAD: `5c925b90e569ae38ef881c0850e0a9be33288b90`
- Upstream: `origin/master`
- Working tree before audit: clean
- Project-level `AGENTS.md` before audit: absent
- Main language: MATLAB
- License: MIT
- MATLAB executable found: `/home/zhuzibn/.local/bin/matlab`
- Source line endings: existing root MATLAB files use CRLF
- Documentation and benchmark MATLAB files use LF

## Purpose and Domain Model

The README describes ASTRO as an atomistic spin-dynamics model developed by
the COSMAL group for 2D GdFeCo ferrimagnets, with possible adaptation to
ferromagnets, antiferromagnets, and 3D structures.

Current code conventions:

- `atomtype_ == 1`: rare-earth Gd site.
- `atomtype_ == 0`: transition-metal FeCo site.
- Spin arrays are indexed `(W, L, time)`.
- Matrix dimension 1 is physical width `W`.
- Matrix dimension 2 is physical length `L`.
- Plotting maps `L` to x and `W` to y, reversing displayed width order.
- `bc == 1` means non-periodic boundaries; `bc == 0` means periodic.
- Default lattice: `20 x 30`, lattice constant `d = 0.4e-9 m`.
- Default composition: `compositionn = 0.1`, meaning 10% Gd sites.

The presentation describes a `0.5 nm` Co0.7Gd0.3 example and a historical
`3 x 5` plotting case. Those values are not the current simulation defaults;
the benchmark correctly preserves the code's `0.4 nm`, 10% Gd behavior.

## Architecture

### Simulation Entry Point

`main.m` is the only complete top-level simulation entry point. It defines
all run parameters in the shared MATLAB workspace, initializes the system,
runs the integrator, and saves the complete workspace to `final.mat`.

This design makes the current behavior easy to run but tightly couples
configuration, state, solver implementation, and output schema.

### Initialization

- `constantfile.m` defines physical constants.
- `systemgeneration.m` calls `distrib.m` and constructs the initial spin
  state, including the optional domain wall or loaded restart state.
- `distrib.m` creates a random RE/TM site distribution with `randperm`, or
  loads/saves the hard-coded `atomtypee.mat`.

Ordinary `main.m` runs call `rng('shuffle')`, so their atom distributions are
not reproducible unless fixed-distribution loading is enabled.

### Effective Field

`field_calc.m` computes:

- nearest-neighbor exchange;
- easy-axis anisotropy;
- interfacial DMI;
- optional thermal field;
- optional dipole field through `dipole_calc.m`; and
- external field.

The default benchmark path has thermal and dipole fields disabled.

### LLG Integration

`integrate_llg.m` allocates GPU arrays, constructs per-site material and
exchange parameters, chunks the simulation, calculates RK4 stages, normalizes
each updated spin, gathers downsampled states, and returns `mmx`, `mmy`,
`mmz`, and `t` in the shared workspace.

Solver kernels:

- `atomgpurk4.m`: one LLG derivative evaluation including SOT/STT terms.
- `atomgpu.m`: Heun update using a fixed supplied effective field.
- `atomgpupc4.m`: fourth-order predictor-corrector update using a fixed
  supplied effective field.

### Dipole Implementations

`dipole_calc.m` has four modes:

- `0`: disabled.
- `1`: direct CPU calculation.
- `2`: direct GPU calculation, explicitly marked potentially wrong.
- `3`: CPU macrocell calculation.

Only mode 0 is covered by the stored regression baseline.

### Post-processing

`plot_dat_proc/` is excluded by `.gitignore` and is not tracked in the current
commit, although files are present in the working directory.

- `threedplott.m` implements the W/L display convention but requires the
  missing external `arrow3` function.
- `threedplott_main.m` and the active section of `data_proc.m` use
  machine-specific paths and expected local MAT-files.

The benchmark figures are self-contained and do not depend on this directory.

## Build, Test, Lint, and Benchmark Commands

Run commands from the repository root.

### Build

No compilation or build command is required or discovered. MATLAB executes
the `.m` files directly.

### Main Simulation

```bash
matlab -batch "run('main.m')"
```

Prerequisites: MATLAB with Parallel Computing Toolbox GPU support, a
compatible NVIDIA GPU, and a writable repository root. This command creates
or overwrites `final.mat` and uses a shuffled atom distribution.

### Quick Regression Benchmark

```bash
matlab -batch "run('benchmark/run_baseline_benchmark.m')"
```

This runs the deterministic `3 x 5` case. It overwrites the shared benchmark
result, figures, and metadata with quick-mode output.

### Current-size Regression Benchmark

```bash
matlab -batch "benchmark_mode='current'; run('benchmark/run_baseline_benchmark.m')"
```

This runs the deterministic `20 x 30` case and overwrites the same benchmark
artifacts with current-mode output.

### Compare a Future Result

```bash
matlab -batch "addpath('benchmark'); passed=compare_to_baseline('path/to/future_result.mat'); assert(passed)"
```

The future result must use the benchmark output schema. The comparison checks
sizes, exact atom distribution, final average spin, final spin arrays, and
spin-norm statistics. It reports but does not include the full-trajectory
maximum difference in the pass condition.

### Tests and Lint

- Automated unit/integration test suite: none discovered.
- MATLAB Code Analyzer command: none documented.
- CI command/configuration: none discovered.
- Octave compatibility: not documented and unlikely because the solver uses
  MATLAB GPU APIs.

## Existing Baseline

The tracked baseline metadata records:

- Date/time: `2026-06-09 17:44:26 +0800`
- Mode: `current`
- Source commit: `9d00a6d735c72a13f23e4477ce09d550e90f1118`
- MATLAB: R2024b Update 1
- GPU: NVIDIA GeForce RTX 3060 Laptop GPU
- Lattice: `20 x 30`
- Saved trajectory: `20 x 30 x 100`
- Runtime: `2e-12 s`
- Solver elapsed time: `149.875465 s`
- Final average spin: approximately `[0.3102, 0.3761, -0.0825]`
- Final maximum spin-norm error: approximately `2.2204e-16`

The audit loaded the MAT-file metadata and selected scalar results without
rerunning or modifying the baseline.

## Risks and Regression Points

1. **Non-default solvers appear broken.** In `integrate_llg.m`, the Heun and
   predictor-corrector branches use `mmxtmp`, `mmytmp`, `mmztmp`, and field
   values before the branch initializes or refreshes them. The
   predictor-corrector bootstrap also calls `atomgpurk4.m` with an obsolete
   argument list and treats derivative outputs as updated spins.
2. **Time/state alignment is unclear.** Each chunk stores `gpusteps` states
   but performs `gpusteps - 1` updates. The next chunk repeats the carried
   boundary state, while `t` is independently generated with `linspace`.
3. **Default runs are non-reproducible.** `main.m` uses `rng('shuffle')` and
   saves a full workspace rather than an explicit result schema.
4. **Workspace-script coupling is high.** Most core modules are scripts that
   depend on many implicit variables. Small parameter or naming changes can
   cause distant failures.
5. **The benchmark duplicates the `main.m` parameter block.** A source change
   can leave benchmark settings stale unless both locations are reviewed.
6. **Thermal-field semantics need validation.** RK4 calls `field_calc.m` four
   times per integration step, so enabling thermal noise generates a new
   random field at each RK4 stage, despite the comment saying once per step.
7. **Dipole modes lack regression coverage.** Mode 2 is explicitly suspect;
   modes 1 and 3 are expensive and not represented in the baseline.
8. **No physics observable for energy exists.** Regression checks are based
   on trajectories, averages, and spin normalization, not energy behavior.
9. **Benchmark pass criteria omit full-trajectory tolerance.** The full
   trajectory difference is printed but only the final spin-array difference
   contributes to `passed`.
10. **Post-processing is not portable.** It contains hard-coded Windows paths,
    expected local result files, inconsistent historical variable names, and
    a missing external plotting dependency.
11. **Source comments are partly stale.** For example, `atomgpurk4.m` is
    labeled as a Heun solver although it returns an RK4 derivative.
12. **No automated branch coverage exists.** Periodic boundaries, thermal
    noise, fixed edges, restart loading, dipole modes, STT, and non-RK4
    solvers are untested.

## Candidate Improvement Areas

Ordered by expected risk reduction:

1. Define and test exact saved-state timestamps and chunk-boundary behavior.
2. Add focused smoke tests for initialization, spin normalization, boundary
   conditions, and field terms using deterministic inputs.
3. Repair or explicitly disable unsupported Heun and predictor-corrector
   paths.
4. Extract a shared parameter/configuration source used by both `main.m` and
   the benchmark without changing default physics.
5. Convert workspace scripts to functions incrementally, starting with pure
   initialization and field calculations.
6. Expand regression coverage to selected disabled terms, especially periodic
   boundaries and deterministic field-only cases.
7. Define an explicit production output schema instead of saving the complete
   workspace.
8. Make plotting portable and self-contained while preserving W/L orientation.
9. Add CI only after a CPU-independent static check or accessible GPU runner
   strategy is defined.

## Files Inspected

- Root: `.gitattributes`, `.gitignore`, `LICENSE`, `README.md`,
  `Readme_improvements.md`
- Simulation: `main.m`, `constantfile.m`, `systemgeneration.m`, `distrib.m`,
  `integrate_llg.m`, `field_calc.m`, `dipole_calc.m`
- Kernels: `atomgpu.m`, `atomgpurk4.m`, `atomgpupc4.m`
- Benchmark: `benchmark/README.md`,
  `benchmark/run_baseline_benchmark.m`,
  `benchmark/compare_to_baseline.m`,
  `benchmark/baseline/baseline_metadata.md`, and the variable schema/scalar
  metadata in `benchmark/baseline/output/baseline_result.mat`
- Documentation: `docs/benchmark_audit.md`,
  `docs/atomistic_model-2D.pptx`
- Present but ignored post-processing:
  `plot_dat_proc/data_proc.m`, `plot_dat_proc/threedplott.m`,
  `plot_dat_proc/threedplott_main.m`
- Git history: recent commits, branches, tracked-file list, and current HEAD

## Recommended Next Phase

Phase 2 should verify the checked-in deterministic baseline on a new branch
before any source improvement. Use the quick mode first because both modes
overwrite the same tracked benchmark artifacts.

Recommended branch:

```text
codex-c3a1-verify-current-baseline
```

Exact next command:

```bash
git switch -c codex-c3a1-verify-current-baseline
```

Then start Codex with:

```bash
codex --sandbox workspace-write --ask-for-approval on-request
```

Exact prompt for the next session:

```text
Phase 2: establish and verify the ASTRO regression baseline.

Read AGENTS.md, audit.md, benchmark/README.md, and
benchmark/baseline/baseline_metadata.md first.

Scope:
- Do not modify MATLAB source physics.
- Focus on benchmark/ and the stored baseline artifacts.
- Run the quick benchmark first:
  matlab -batch "run('benchmark/run_baseline_benchmark.m')"
- Record the exact environment, command, runtime, output schema, and result.
- Compare the generated result with the accepted current-size baseline only
  if a separate future-result file can be produced without overwriting the
  accepted artifact.
- Do not install dependencies, use the network, commit, or push.
- Stop before committing.
- Report commands, results, failures, and whether the baseline should be
  committed.
```
