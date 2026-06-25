# ASTRO Project Instructions

## Project Purpose

ASTRO is a MATLAB atomistic spin-dynamics model for a 2D ferrimagnetic
Gd-FeCo lattice. The current production path runs on an NVIDIA GPU through
MATLAB GPU arrays and integrates the Landau-Lifshitz-Gilbert equation.

Read `audit.md` for the current architecture, risks, and audit evidence.

## Important Entry Points

- `main.m`: top-level simulation and parameter definition.
- `constantfile.m`: physical constants.
- `systemgeneration.m`: atom distribution and initial spin state.
- `distrib.m`: random or fixed RE/TM site distribution.
- `integrate_llg.m`: GPU allocation, exchange setup, time integration, and
  saved trajectory assembly.
- `field_calc.m`: exchange, anisotropy, DMI, thermal, dipole, and external
  effective fields.
- `dipole_calc.m`: dipole modes 0-3.
- `atomgpurk4.m`: LLG derivative used by the active RK4 path.
- `atomgpu.m`: Heun kernel.
- `atomgpupc4.m`: predictor-corrector kernel.
- `benchmark/run_baseline_benchmark.m`: deterministic regression runner.
- `benchmark/compare_to_baseline.m`: benchmark result comparison.

`plot_dat_proc/` is ignored and contains non-portable historical
post-processing. Do not treat it as maintained production code without
explicit approval.

## Environment

- Run commands from the repository root.
- MATLAB with Parallel Computing Toolbox GPU support is required.
- A compatible NVIDIA GPU is required for the current solver path.
- `rk4=1` is the only supported production solver mode. `rk4=0` and
  `rk4=2` must fail with `ASTRO:UnsupportedSolver`; do not repair or use the
  Heun or predictor-corrector paths without explicit validation scope.
- Do not assume GNU Octave compatibility.
- Do not add dependencies or use the network without approval.

## Commands

Main simulation:

```bash
matlab -batch "run('main.m')"
```

This uses `rng('shuffle')` and creates or overwrites root `final.mat`.

Quick deterministic benchmark:

```bash
matlab -batch "run('benchmark/run_baseline_benchmark.m')"
```

Current-size deterministic benchmark:

```bash
matlab -batch "benchmark_mode='current'; run('benchmark/run_baseline_benchmark.m')"
```

Compare a future benchmark result:

```bash
matlab -batch "addpath('benchmark'); passed=compare_to_baseline('path/to/future_result.mat'); assert(passed)"
```

Deterministic invariant smoke checks:

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
```

No build command, unit-test framework, lint command, or CI workflow is
currently defined.

## Numerical and Domain Conventions

- `atomtype_ == 1` means RE/Gd.
- `atomtype_ == 0` means TM/FeCo.
- Arrays are indexed `(W, L, time)`.
- First matrix dimension is physical width `W`.
- Second matrix dimension is physical length `L`.
- Plot physical length on x and physical width on y.
- Reverse displayed width order to preserve the positive physical y
  direction used by the model.
- `bc == 1` is non-periodic; `bc == 0` is periodic.
- Current code uses `d = 0.4e-9 m` and `compositionn = 0.1`.
- Do not replace current parameters with values from the presentation unless
  the task explicitly requests a physics change.
- Spin normalization is an important invariant after each integration update.
- GPU results may require tolerances across MATLAB, CUDA, driver, and GPU
  versions; do not require bitwise identity without evidence.

## Benchmark and Baseline Rules

- Treat `benchmark/baseline/` as accepted regression data.
- Read `benchmark/README.md` and baseline metadata before changing benchmark
  behavior.
- Both quick and current benchmark modes overwrite the same result, figures,
  and metadata. Preserve accepted artifacts before exploratory runs.
- Do not update baseline artifacts as a side effect of unrelated work.
- Do not change comparison tolerances without recorded cross-platform
  evidence and explicit approval.
- The baseline protects current behavior; it is not proof of physical
  correctness.
- No energy observable currently exists. Do not invent one in benchmark
  documentation or results.
- Keep fixed atom distributions versioned and preserve the meaning
  `1=RE/Gd`, `0=TM/FeCo`.

## Files Requiring Approval

Do not modify these without explicit task scope:

- `benchmark/baseline/input/*.mat`
- `benchmark/baseline/output/*.mat`
- `benchmark/baseline/figures/*`
- `benchmark/baseline/baseline_metadata.md`
- `docs/atomistic_model-2D.pptx`
- `LICENSE`

Do not rewrite Git history, delete result data, or overwrite root `final.mat`
unless the user approves the operation.

## Code Style

- Preserve existing line endings. Root MATLAB source currently uses CRLF.
- Use small patch-based edits and avoid formatting-only churn.
- Match existing MATLAB naming and script style unless a scoped refactor
  explicitly changes it.
- Keep units in comments for physical parameters.
- State whether a value is per atom, per link, field, moment, length, or time
  when adding parameters.
- Avoid new workspace dependencies. Prefer explicit function inputs and
  outputs for new isolated code.
- Do not silently change physics formulas, signs, axis orientation, boundary
  semantics, solver selection, time indexing, random distributions, or
  normalization.
- Preserve explicit rejection of unsupported solver modes unless a task
  explicitly validates and restores those solver paths.

## Documentation

- Update root `CHANGELOG.md` for source features, bug fixes, or resolved
  non-trivial errors, following the global instructions.
- Do not update the changelog for audit-only, explanation-only, or
  documentation-only work.
- Update `audit.md` or project instructions only when repository facts change.
- Mark uncertain commands or behavior as tentative rather than guessing.

## Repository Modification Workflow

1. Audit:
   - Check `git status --short --branch`.
   - Stop and ask if the working tree is dirty.
   - Read `AGENTS.md`, `audit.md`, relevant source, benchmark docs, and recent
     history.
   - Identify the exact behavior, inputs, outputs, and regression risks.
2. Baseline:
   - Choose quick or current benchmark deliberately.
   - Record the exact command, MATLAB/GPU environment, runtime, and artifacts.
   - Avoid overwriting the accepted baseline unless that is the task.
3. Plan:
   - Split work into one behavioral aspect per step.
   - State files, risk, expected benefit, verification, branch, and commit
     recommendation for each step.
4. Branch:
   - Use `codex-<short-id>-<verb>-<topic>`.
   - Create a branch before edits unless the user explicitly requests direct
     edits.
5. Implement:
   - Make the smallest change that satisfies the requested behavior.
   - Do not refactor adjacent code or alter disabled physics paths unless in
     scope.
6. Verify:
   - Run the narrowest relevant deterministic benchmark or check.
   - Check spin shape, finite values, normalization, atom distribution, and
     saved timestamps when relevant.
   - Run `git diff --stat` and `git diff --check`.
   - Investigate unexpected large or formatting-only diffs.
7. Report:
   - List branch, changed files, commands, test result, benchmark result,
     commit readiness, suggested commit command, and next step.
   - Do not commit or push unless explicitly requested.

## Known High-Risk Areas

- Heun and predictor-corrector branches are not covered by the baseline and
  appear inconsistent with current variable initialization/function
  signatures.
- Saved-state timestamps and chunk-boundary duplication need explicit
  validation before changes.
- Thermal noise is generated during every RK4 field evaluation when enabled.
- Dipole mode 2 is marked potentially wrong in source comments.
- The benchmark duplicates parameters from `main.m`; review both when default
  physics parameters change.
- Post-processing uses hard-coded paths and a missing external `arrow3`
  dependency.
