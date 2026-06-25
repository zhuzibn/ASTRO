# ASTRO Step-by-Step Improvement Plan

Plan date: 2026-06-10

Phase: 3 - improve step by step.

This plan is based on `audit.md`, the accepted current-size regression
baseline, and the current `master` branch at `1013a9f`. It does not authorize
source edits, baseline replacement, commits, or pushes by itself.

## Current Preconditions

- The accepted baseline is `benchmark/baseline/output/baseline_result.mat`.
- The accepted baseline was produced on MATLAB R2024b Update 1 with an NVIDIA
  GeForce RTX 3060 Laptop GPU.
- There is no unit-test framework, lint command, CI workflow, or CPU solver.
- The current working tree contains a pre-existing deletion of
  `docs/benchmark_audit.md`. Do not include, restore, or commit that deletion
  as part of any step unless separately instructed.
- Root MATLAB source uses CRLF. Preserve existing line endings.
- Every source-code fix or feature must update root `CHANGELOG.md`. Since it
  does not currently exist, the first source-changing step must initialize it
  from `~/.config/opencode/CHANGELOG.md`, or create a minimal equivalent if
  the template is unavailable.
- Before starting an implementation branch, obtain a clean accepted starting
  point or explicit permission to carry unrelated working-tree changes.
- Do not modify accepted baseline inputs, output, figures, metadata,
  `docs/atomistic_model-2D.pptx`, or `LICENSE` without explicit approval.

## Sequence and Decision Gates

Execute the steps in order. Each step addresses one primary concern and
should be reviewed before the next begins.

1. Make exploratory benchmark runs non-destructive.
2. Define and correct the saved-state/time contract.
3. Add deterministic invariant smoke checks.
4. Fail clearly for unsupported solver selections.
5. Remove duplicated production/benchmark parameter definitions.
6. Isolate deterministic initialization from workspace scripts.
7. Add deterministic boundary and field-term coverage.
8. Define an explicit production result schema.

Do not repair Heun or predictor-corrector numerics until their intended field
refresh, torque support, startup method, and validation criteria are defined.
Step 4 deliberately disables unsupported selections instead of guessing at
their physics.

## Step 1 - Non-Destructive Benchmark Output

**Goal:** Allow quick and current benchmark runs to write candidate results,
figures, and metadata outside `benchmark/baseline/`, while preserving the
current default baseline-generation behavior only when explicitly requested.

**Likely files:** `benchmark/run_baseline_benchmark.m`,
`benchmark/README.md`, `CHANGELOG.md`.

**Risk:** Low. The solver and physics parameters should remain unchanged.

**Expected benefit:** Later steps can generate comparable results without
overwriting accepted artifacts.

**Verification:**

```bash
matlab -batch "benchmark_output_root=fullfile(tempdir,'astro-step1-quick'); run('benchmark/run_baseline_benchmark.m')"
matlab -batch "addpath('benchmark'); passed=compare_to_baseline(fullfile(tempdir,'astro-step1-quick','output','benchmark_result.mat')); assert(passed)"
git diff --stat
git diff --check
```

The exact candidate filename may be adjusted during implementation, but it
must be documented and the accepted baseline tree must remain byte-for-byte
unchanged.

**Branch:** New branch from the accepted `master` baseline.

```bash
git status --short --branch
git switch -c codex-d1a1-add-candidate-benchmark-output
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit separately after the quick candidate run
passes comparison and accepted baseline artifacts are unchanged.

```bash
git add benchmark/run_baseline_benchmark.m benchmark/README.md CHANGELOG.md
git commit -m "Add non-destructive benchmark output mode"
```

**Codex prompt:**

```text
Phase 3 implementation step 1: add non-destructive ASTRO benchmark output.

Read AGENTS.md, audit.md, fix.md, benchmark/README.md, and the benchmark
runner first. Preserve the accepted benchmark/baseline/ tree exactly.

Scope:
- Add one documented way to direct benchmark result, figures, and metadata
  to a caller-selected candidate directory.
- Do not change solver behavior, physics parameters, comparison tolerances,
  fixed atom distributions, or accepted baseline artifacts.
- Initialize/update root CHANGELOG.md as required by AGENTS.md.
- Run a quick candidate benchmark in tempdir and compare it to the accepted
  baseline when dimensions/mode are compatible.
- Run git diff --stat and git diff --check.
- Do not commit or push.
- Report exact commands, outputs, changed files, and commit readiness.
```

## Step 2 - Saved-State and Time Contract

**Goal:** Define `mmx(:,:,k)`, `mmy(:,:,k)`, `mmz(:,:,k)`, and `t(k)` as the
same physical state and time, include the initial state at `t=0`, include the
final state at `runtime`, and avoid duplicate chunk-boundary states.

**Likely files:** `main.m`, `integrate_llg.m`,
`benchmark/run_baseline_benchmark.m`, `benchmark/compare_to_baseline.m`,
`benchmark/README.md`, `CHANGELOG.md`.

**Risk:** High. This intentionally changes trajectory shape, timestamps, and
accepted behavior. Do not replace the accepted baseline in this step.

**Expected benefit:** Removes the largest known ambiguity in simulation
output and makes trajectory comparisons physically interpretable.

**Verification:**

```bash
matlab -batch "benchmark_output_root=fullfile(tempdir,'astro-step2-quick'); run('benchmark/run_baseline_benchmark.m')"
matlab -batch "r=load(fullfile(tempdir,'astro-step2-quick','output','benchmark_result.mat')); assert(numel(r.t)==size(r.mmx,3)); assert(r.t(1)==0); assert(abs(r.t(end)-r.runtime)<=eps(r.runtime)); assert(all(diff(r.t)>0)); assert(all(isfinite(r.mmx(:)))); assert(all(isfinite(r.mmy(:)))); assert(all(isfinite(r.mmz(:))))"
git diff --stat
git diff --check
```

Also record the expected saved count from `runtime`, `tstep`, and
`savetstep`, and explicitly verify that no chunk boundary appears twice.

**Branch:** New branch from accepted `master`, or a continuation from the
accepted Step 1 commit.

```bash
git status --short --branch
git switch -c codex-d1a2-fix-saved-time-grid
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit after the new contract passes the quick
checks. Baseline regeneration must be a separate approved step.

```bash
git add main.m integrate_llg.m benchmark/run_baseline_benchmark.m benchmark/compare_to_baseline.m benchmark/README.md CHANGELOG.md
git commit -m "Align saved spin states with simulation times"
```

**Codex prompt:**

```text
Phase 3 implementation step 2: fix ASTRO saved-state/time alignment.

Read AGENTS.md, audit.md, fix.md, main.m, integrate_llg.m, and benchmark
documentation first. Use the non-destructive benchmark output mechanism from
step 1.

Required contract:
- The first saved state is the initial state at t=0.
- The last saved state is the state at runtime.
- numel(t) equals size(mmx,3), size(mmy,3), and size(mmz,3).
- Saved times are strictly increasing.
- Chunk boundaries are not duplicated.

Make the smallest RK4-path change that satisfies this contract. Do not alter
field formulas, torque formulas, boundary semantics, normalization, random
distribution behavior, or accepted baseline artifacts. Update CHANGELOG.md.
Run the quick candidate benchmark and explicit assertions for shape, finite
values, normalization, endpoints, and monotonic timestamps. Run git diff
--stat and git diff --check. Do not commit or push. Report exact commands and
results.
```

## Step 3 - Deterministic Invariant Smoke Checks

**Goal:** Add a small MATLAB smoke-check entry point for deterministic
initialization, output shape, finite values, spin normalization, atom counts,
and saved timestamps.

**Likely files:** New `tests/run_smoke_tests.m` or
`benchmark/run_smoke_tests.m`, benchmark helper code if reuse is necessary,
`README.md`, `AGENTS.md`, `CHANGELOG.md`.

**Risk:** Low to medium. Tests must not mutate accepted baseline artifacts.

**Expected benefit:** Provides a repeatable verification command narrower
than the current-size benchmark and catches structural regressions early.

**Verification:**

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
git diff --stat
git diff --check
```

Use the discovered path if the implementation chooses `benchmark/` instead
of `tests/`. The test must exit nonzero on failure.

**Branch:** New branch from the accepted Step 2 result.

```bash
git status --short --branch
git switch -c codex-d1a3-add-invariant-smoke-tests
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit separately when the smoke command passes.

```bash
git add tests README.md AGENTS.md CHANGELOG.md
git commit -m "Add deterministic ASTRO smoke checks"
```

Adjust the path list if tests are placed under `benchmark/`.

**Codex prompt:**

```text
Phase 3 implementation step 3: add deterministic ASTRO invariant smoke
checks.

Read AGENTS.md, audit.md, fix.md, benchmark/README.md, and the accepted
saved-state/time contract first.

Scope:
- Add one documented MATLAB batch command that exits nonzero on failure.
- Check deterministic atom distribution/counts, trajectory dimensions,
  finite values, spin normalization, t=0, t(end)=runtime, strictly increasing
  saved times, and matching time/state counts.
- Reuse existing benchmark inputs and helpers where practical, without
  writing to benchmark/baseline/.
- Do not add dependencies or change physics.
- Update README.md, AGENTS.md commands, and CHANGELOG.md.
- Run the smoke command, git diff --stat, and git diff --check.
- Do not commit or push. Report exact commands and results.
```

## Step 4 - Explicitly Reject Unsupported Solvers

**Goal:** Reject `rk4=0` and `rk4=2` before GPU integration with a clear
message explaining that only the RK4 path is currently supported.

**Likely files:** `main.m`, `integrate_llg.m`,
`benchmark/run_baseline_benchmark.m`, smoke tests, `README.md`,
`CHANGELOG.md`. Do not modify `atomgpu.m` or `atomgpupc4.m` in this step.

**Risk:** Low. Previously broken selections become explicit errors.

**Expected benefit:** Prevents undefined-variable and obsolete-signature
failures from being mistaken for supported solver behavior.

**Verification:**

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
matlab -batch "rk4=0; try, error('Test harness must invoke solver validation here'); catch e, disp(e.message); end"
git diff --stat
git diff --check
```

The implementation must replace the illustrative second command with a real
test that asserts both unsupported values fail for the intended reason.

**Branch:** New branch from the accepted Step 3 result.

```bash
git status --short --branch
git switch -c codex-d1a4-disable-unsupported-solvers
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit separately after positive RK4 and negative
solver-selection checks pass.

```bash
git add main.m integrate_llg.m benchmark/run_baseline_benchmark.m tests README.md CHANGELOG.md
git commit -m "Reject unsupported ASTRO solver modes"
```

**Codex prompt:**

```text
Phase 3 implementation step 4: explicitly reject unsupported ASTRO solver
modes.

Read AGENTS.md, audit.md, fix.md, main.m, integrate_llg.m, atomgpu.m,
atomgpupc4.m, and current smoke tests.

Scope:
- Treat rk4=1 as the only supported production solver.
- Make rk4=0 and rk4=2 fail before integration with a precise message.
- Add deterministic negative checks for both unsupported values.
- Do not attempt to repair Heun or predictor-corrector and do not edit their
  kernels.
- Do not change RK4 numerics or accepted baseline artifacts.
- Update relevant documentation and CHANGELOG.md.
- Run smoke tests, git diff --stat, and git diff --check.
- Do not commit or push. Report exact commands and results.
```

## Step 5 - Shared Default Configuration

**Goal:** Define production physics and time defaults once, then consume the
same values from `main.m` and the benchmark while allowing only documented
benchmark overrides such as lattice size and fixed atom distribution.

**Likely files:** New `astro_default_config.m`, `main.m`,
`benchmark/run_baseline_benchmark.m`, smoke tests, `README.md`,
`AGENTS.md`, `CHANGELOG.md`.

**Risk:** Medium to high. A missed or renamed workspace variable can silently
change behavior.

**Expected benefit:** Prevents drift between production defaults and the
regression runner.

**Verification:**

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
matlab -batch "benchmark_output_root=fullfile(tempdir,'astro-step5-quick'); run('benchmark/run_baseline_benchmark.m')"
matlab -batch "addpath('benchmark'); passed=compare_to_baseline(fullfile(tempdir,'astro-step5-quick','output','benchmark_result.mat')); assert(passed)"
git diff --stat
git diff --check
```

Comparison should be against the post-time-contract accepted candidate, not
an incompatible pre-contract trajectory.

**Branch:** New branch from the accepted Step 4 result.

```bash
git status --short --branch
git switch -c codex-d1a5-share-default-configuration
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit only if parameter snapshots and trajectory
comparison show no unintended behavior change.

```bash
git add astro_default_config.m main.m benchmark/run_baseline_benchmark.m tests README.md AGENTS.md CHANGELOG.md
git commit -m "Share ASTRO default simulation configuration"
```

**Codex prompt:**

```text
Phase 3 implementation step 5: remove duplicated ASTRO default parameter
definitions.

Read AGENTS.md, audit.md, fix.md, main.m, the benchmark runner, and current
tests first.

Scope:
- Add one simple explicit default-configuration function or struct.
- Make main.m and the benchmark consume the same production defaults.
- Keep benchmark-only overrides limited to documented lattice size, fixed
  atom distribution, and candidate output paths.
- Preserve every current default value, unit, enabled term, and solver
  selection.
- Avoid a general framework or unrelated script-to-function conversion.
- Update tests, README.md, AGENTS.md if commands change, and CHANGELOG.md.
- Run smoke tests and a non-destructive trajectory comparison.
- Run git diff --stat and git diff --check.
- Do not commit or push. Report exact commands and results.
- suggest do I need to merge the changes to the master branch? do I need to push to github?
```

## Step 6 - Pure Deterministic Initialization

**Goal:** Extract atom-distribution validation and initial-spin construction
into explicit functions with inputs and outputs, while leaving random/fixed
distribution policy in the caller.

**Likely files:** New initialization function files, `systemgeneration.m`,
`benchmark/run_baseline_benchmark.m`, smoke tests, `README.md`,
`CHANGELOG.md`.

**Risk:** Medium. RE/TM sign conventions and domain-wall orientation must not
change.

**Expected benefit:** Makes initialization independently testable and reduces
workspace coupling without touching integration or field physics.

**Verification:**

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
matlab -batch "benchmark_output_root=fullfile(tempdir,'astro-step6-quick'); run('benchmark/run_baseline_benchmark.m')"
git diff --stat
git diff --check
```

Tests must compare extracted initialization output exactly against the
previous deterministic fixed-input result.

**Branch:** New branch from the accepted Step 5 result.

```bash
git status --short --branch
git switch -c codex-d1a6-extract-spin-initialization
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit separately after exact initialization and
trajectory checks pass.

```bash
git add systemgeneration.m benchmark/run_baseline_benchmark.m tests README.md CHANGELOG.md
git add <new-initialization-function-files>
git commit -m "Extract deterministic spin initialization"
```

**Codex prompt:**

```text
Phase 3 implementation step 6: extract pure deterministic ASTRO
initialization.

Read AGENTS.md, audit.md, fix.md, systemgeneration.m, distrib.m, the benchmark
runner, and current tests.

Scope:
- Extract only fixed atom-distribution validation and initial spin
  construction into explicit functions.
- Preserve atomtype_ meaning, W/L indexing, domain-wall orientation, angles,
  and exact deterministic output.
- Keep random/load/save policy in existing callers.
- Do not modify integration, fields, solver kernels, or physics defaults.
- Add exact-output tests, update documentation and CHANGELOG.md.
- Run smoke tests, a quick candidate benchmark, git diff --stat, and git diff
  --check.
- Do not commit or push. Report exact commands and results.
```

## Step 7 - Boundary and Deterministic Field Coverage

**Goal:** Add focused deterministic checks for non-periodic and periodic
neighbor semantics plus exchange, anisotropy, DMI, and external-field terms
on very small hand-verifiable lattices.

**Likely files:** New test/helper files, possibly a narrowly extracted
exchange-setup or field function, `field_calc.m`, `integrate_llg.m`,
`README.md`, `AGENTS.md`, `CHANGELOG.md`.

**Risk:** High. Field signs, dimensions, and boundary orientation are
physics-sensitive.

**Expected benefit:** Adds coverage for core formulas that the trajectory
baseline alone cannot localize.

**Verification:**

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
matlab -batch "run('tests/run_field_tests.m')"
git diff --stat
git diff --check
```

Expected values must be derived and documented from the existing formulas.
Do not add thermal or dipole coverage in this step.

**Branch:** New branch from the accepted Step 6 result.

```bash
git status --short --branch
git switch -c codex-d1a7-test-boundaries-and-fields
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit separately only after hand-derived expected
values and both boundary modes pass.

```bash
git add tests field_calc.m integrate_llg.m README.md AGENTS.md CHANGELOG.md
git add <new-field-helper-files>
git commit -m "Add deterministic boundary and field checks"
```

**Codex prompt:**

```text
Phase 3 implementation step 7: add deterministic boundary and field-term
coverage.

Read AGENTS.md, audit.md, fix.md, integrate_llg.m, field_calc.m,
constantfile.m, and current tests first.

Scope:
- Test bc=1 and bc=0 neighbor behavior on small fixed atom matrices.
- Test exchange, anisotropy, DMI, and external field separately with
  hand-derived expected values.
- Extract only the minimum pure helper needed for testing.
- Preserve formulas, signs, W/L orientation, units, and GPU production path.
- Do not enable or modify thermal fields, dipole modes, torque formulas, or
  accepted baseline artifacts.
- Document derivations succinctly and update CHANGELOG.md.
- Run smoke and field tests, git diff --stat, and git diff --check.
- Do not commit or push. Report exact commands and results.
```

## Step 8 - Explicit Production Result Schema

**Goal:** Replace `save('final.mat')` of the full workspace with an explicit
versioned result struct or named-variable schema, while retaining an explicit
compatibility option only if required.

**Likely files:** `main.m`, new result-save helper if useful, tests,
`README.md`, `AGENTS.md`, `CHANGELOG.md`.

**Risk:** Medium to high. Existing downstream scripts may rely on undeclared
workspace variables.

**Expected benefit:** Makes production outputs stable, documented, smaller,
and auditable.

**Verification:**

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
matlab -batch "result_file=fullfile(tempdir,'astro-final.mat'); run('main.m'); r=whos('-file',result_file); assert(~isempty(r))"
git diff --stat
git diff --check
```

The actual main-run command must avoid overwriting root `final.mat` and may
require a documented output-path input added by this step.

**Branch:** New branch from the accepted Step 7 result.

```bash
git status --short --branch
git switch -c codex-d1a8-define-result-schema
codex --sandbox workspace-write --ask-for-approval on-request
```

**Commit recommendation:** Commit separately after schema assertions and a
non-destructive main run pass. Do not delete historical result files.

```bash
git add main.m tests README.md AGENTS.md CHANGELOG.md
git add <new-result-helper-files>
git commit -m "Define explicit ASTRO result schema"
```

**Codex prompt:**

```text
Phase 3 implementation step 8: define an explicit ASTRO production result
schema.

Read AGENTS.md, audit.md, fix.md, main.m, benchmark output schema, and current
tests first.

Scope:
- Add a caller-selectable production output path.
- Save an explicit documented versioned schema instead of the complete MATLAB
  workspace.
- Include configuration, atom distribution, initial state, trajectory,
  timestamps, and essential provenance.
- Do not overwrite root final.mat during verification.
- Preserve simulation physics and trajectory values.
- Do not delete data or update accepted benchmark artifacts.
- Add schema tests and update README.md, AGENTS.md, and CHANGELOG.md.
- Run smoke tests, a non-destructive main run, git diff --stat, and git diff
  --check.
- Do not commit or push. Report exact commands and results.
```

## Deferred Work

The following items need separate design and validation plans after the eight
steps above:

- Repairing Heun or fourth-order predictor-corrector integration.
- Defining thermal-noise sampling semantics for RK4 stages.
- Validating dipole modes 1-3, especially mode 2.
- Adding energy observables.
- Making ignored historical post-processing portable.
- Adding CI before a suitable MATLAB/GPU runner or meaningful static-only
  check is available.

## Phase 3 Completion Criteria

- `fix.md` is reviewed and accepted.
- No source, physics, benchmark baseline, or global instruction file changed
  during planning.
- Every implementation step has one primary concern, files, risk, benefit,
  verification, branch guidance, commit guidance, commands, and a paste-ready
  Codex prompt.
- Implementation begins only after the working-tree state and starting branch
  for Step 1 are explicitly accepted.
