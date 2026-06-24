# ASTRO Baseline Metadata

- Benchmark date/time: `2026-06-24 16:58:03 +0800`
- Git commit: `638325cfc905ee57ce98f17ae57e6d43cf44a9f1`
- Git branch: `master`
- MATLAB version: `24.2.0.2740171 (R2024b) Update 1`
- GPU available: `1`
- GPU: `NVIDIA GeForce RTX 3060 Laptop GPU`
- GPU compute capability: `8.6`
- GPU driver: `581.95`
- CUDA toolkit version reported by MATLAB: `12.2000`
- Operating system: `Microsoft Windows 11 家庭版 中文版`
- Benchmark mode: `current`
- Exact MATLAB expression: `benchmark_mode='current'; run('benchmark/run_baseline_benchmark.m')`
- Lattice size: `20 x 30` (W x L)
- Lattice constant: `0.4e-9 m`
- Time settings: `tstep=2e-16 s`, `gpusave=1e-12 s`, `gpurun_number=2`, `runtime=2e-12 s`, `savetstep=100`
- Composition: `compositionn=0.1` (10% RE/Gd)
- Fixed atom distribution: `C:\Users\zzf-m\OneDrive\code_softwares\atomistic_model\2D\atomistic\atomistic_v2\benchmark\baseline\input\baseline_atomtype.mat`
- Enabled: exchange, anisotropy, DMI, SOT damping-like torque
- Disabled/zero: thermal, dipole, SOT field-like, both STT terms, fixed edges, external field
- Result file: `C:\Users\zzf-m\OneDrive\code_softwares\atomistic_model\2D\atomistic\atomistic_v2\benchmark\baseline\output\baseline_result.mat`
- Figure directory: `C:\Users\zzf-m\OneDrive\code_softwares\atomistic_model\2D\atomistic\atomistic_v2\benchmark\baseline\figures`
- Benchmark actually ran: `1`
- Solver runtime: `138.972184 seconds`
- Warnings: none recorded
- Missing quantities: energy is not computed by current ASTRO output.

## Deviations From `main.m`

- The runner does not execute `main.m` because it starts with `clear all` and `rng('shuffle')`.
- It copies the audited parameter block, runs the existing `systemgeneration.m`, then replaces the random `atomtype_` with the fixed input and rebuilds the same initial-state formula.
- It calls the unchanged `integrate_llg.m`, `field_calc.m`, and `dipole_calc.m` path.
- It saves an explicit benchmark MAT-file instead of the root-level `final.mat`; no `final.mat` is required or overwritten.
- Quick mode changes only lattice size from `20 x 30` to `3 x 5`; all listed physics and time settings remain unchanged.

## Interpretation

This is a regression baseline for current code behavior. It is not a proof or full validation of the physics.
