# ASTRO Baseline Benchmark Audit

Audit date: 2026-06-09

Scope: repository inspection only. No MATLAB source or benchmark script was
changed or created.

## 1. Current runnable 2D ASTRO entry points

### `main.m`

`main.m` is the only current top-level simulation script found. It:

1. loads constants from `constantfile.m`;
2. creates the atom distribution and initial spins through
   `systemgeneration.m` and `distrib.m`;
3. defines the material, field, torque, and time parameters;
4. runs `integrate_llg.m`; and
5. saves the full workspace as `final.mat`.

The execution path is GPU-based. `main.m` states that an NVIDIA GPU is
required, and `integrate_llg.m` creates `gpuArray` data and uses GPU
`arrayfun`.

### Other scripts

- `plot_dat_proc/threedplott_main.m` is a post-processing entry point, not a
  simulation. It loads a hard-coded external MAT-file path.
- `plot_dat_proc/data_proc.m` is post-processing and comparison code, not a
  clean simulation entry point. Its enabled section expects `final3.mat` and
  `final4.mat`.
- `systemgeneration.m`, `integrate_llg.m`, `field_calc.m`, and
  `dipole_calc.m` are workspace scripts called by `main.m`; they are not
  standalone examples.
- `tmp.m` is an incomplete/debug calculation fragment, not an entry point.

The `plot_dat_proc/` directory and `tmp.m` are present in the inspected
working tree but are currently untracked by Git.

## 2. Best baseline candidate

The safest source of current ASTRO behavior is `main.m`, because it is the
only complete, tracked, top-level 2D simulation path and it exercises the
current `systemgeneration` -> `integrate_llg` -> `field_calc` implementation.

It is not reproducible as currently configured because it executes
`rng('shuffle')` before generating `atomtype_`. A future baseline benchmark
should therefore preserve the `main.m` physics path but explicitly pin the
random seed before `distrib.m`, or load a version-controlled fixed
`atomtype_`. Because `main.m` begins with `clear all`, a seed cannot reliably
be injected by a preceding command without changing or wrapping the current
entry path.

For the first baseline, the lowest-risk configuration is the current
deterministic-physics subset: thermal field disabled and dipole field
disabled, with a fixed atom distribution. The baseline should record both
the parameters and `atomtype_`; recording a seed alone may not guarantee the
same `randperm` result across MATLAB releases.

## 3. Required input parameters

### Geometry and time

| Input | Current value | Meaning |
|---|---:|---|
| `natomW` | `20` | Number of atoms/cells along width; first matrix dimension |
| `natomL` | `30` | Number of atoms/cells along length; second matrix dimension |
| `d` | `0.4e-9 m` | Lattice constant |
| `tstep` | `2e-16 s` | LLG integration timestep |
| `gpusave` | `1e-12 s` | Duration allocated to each GPU run chunk |
| `gpurun_number` | `2` | Number of GPU run chunks |
| `runtime` | `2e-12 s` | Computed as `gpurun_number*gpusave` |
| `savetstep` | `100` | Save every 100 integration storage positions |

Boundary condition is `bc=1`, documented in the code as non-periodic.
The solver is `rk4=1`.

### Random seed and atom distribution

- Current behavior: `rng('shuffle')`.
- Distribution: `randperm(natom, round(natom*compositionn))`.
- Fixed distribution loading: disabled (`load_fixed_atom_distrib=0`).
- Fixed distribution saving: disabled (`save_fixed_atom_distrib=0`).
- Reproducible seed value: **UNCLEAR**; none is defined.
- A fixed distribution file name is `atomtypee.mat`, but no such tracked
  baseline input file was found.

### Co/Gd composition

- `compositionn=0.1`.
- `atomtype_ == 1` means RE/Gd.
- `atomtype_ == 0` means TM/FeCo.
- Therefore the modeled lattice is 10% Gd and 90% FeCo/TM, subject to
  `round(natom*compositionn)`.
- The separate Co fraction within FeCo is **UNCLEAR**; no Fe-versus-Co atom
  distinction or Co concentration parameter was found.

### Material parameters

| Parameter | Current value |
|---|---:|
| Easy-axis anisotropy `Ksim` | `0.4e-3*ele J` |
| Gd-Gd exchange `Jgdgd` | `-1.26e-21 J/link` |
| FeCo-FeCo exchange `Jfefe` | `-2.835e-21 J/link` |
| FeCo-Gd exchange `Jfegd` | `1.09e-21 J/link` |
| TM g-factor `gTM` | `2.2` |
| RE g-factor `gRE` | `2` |
| RE moment `musRE` | `7.63*mub J/T` |
| TM moment `musTM` | `2.217*mub J/T` |
| Film thickness `tz` | `d` |
| Gilbert damping `alp` | `0.02` |
| Temperature `T` | `100 K` (thermal field is disabled) |

Derived gyromagnetic ratios, saturation magnetizations, and torque fields are
also saved in the full-workspace output.

### Enabled fields and torques

| Contribution | Current state | Controlling values |
|---|---|---|
| Nearest-neighbor exchange | Enabled | `Jgdgd`, `Jfefe`, `Jfegd` |
| Uniaxial anisotropy | Enabled | `Ksim` |
| DMI | Enabled | `DMIenable=1`, `Dsim=128e-6*ele J` |
| External field | Numerically zero | `Hext=[0,0,0] T` |
| Thermal field | Disabled | `thermalenable=0` |
| Dipole field | Disabled | `dipolemode=0` |
| SOT damping-like torque | Enabled | `SOT_DLT=1`, `jcSOT=1e9 A/m^2`, `thetaSH=0.2`, `psjSHE=[0,1,0]` |
| SOT field-like torque | Disabled | `SOT_FLT=0`, `chi=0` |
| STT damping-like torque | Disabled | `STT_DLT=0` |
| STT field-like torque | Disabled | `STT_FLT=0`, `chiSTT=0` |
| Fixed edge spins | Disabled | `enablefixedge=0` |

The initial condition uses `dwcalc=1`, which creates a domain wall near the
middle of the length dimension. Loading a prior spin state is disabled with
`loadstartm=0`.

## 4. Current outputs

`main.m` calls `save('final.mat')` without a variable list, so it saves the
entire remaining MATLAB workspace.

### Spin states

- `mx_init`, `my_init`, `mz_init`: initial spin components, each
  `natomW x natomL`.
- `mmx`, `mmy`, `mmz`: saved spin trajectories, each
  `natomW x natomL x Nsaved`.
- `t`: downsampled time vector after `integrate_llg.m`.

For the current settings, the intended saved trajectory size appears to be
`20 x 30 x 100`. Exact time/state alignment should be verified before using
it as a reference; see risks below.

### Atom type matrix

- `atomtype_`: `natomW x natomL`, with `1` for RE/Gd and `0` for TM/FeCo.

### Magnetization

No explicit total or sublattice magnetization variable is produced by
`main.m`.

`plot_dat_proc/data_proc.m` contains examples that derive:

- average TM spin components by summing `mmx`, `mmy`, and `mmz` over TM
  sites; and
- a net `z` magnetization-like quantity named `mnet`, weighted by
  `musRE/mub` or `musTM/mub`.

These are post-processing examples and are not part of the current simulation
output contract. Some disabled `mnet` code uses inconsistent variable names
or indexing, so its suitability as a benchmark observable is **UNCLEAR**.

### Energy

No exchange, anisotropy, DMI, Zeeman, dipole, or total energy output was
found. Energy availability is therefore **UNCLEAR / not currently
implemented as an output**.

### Figures

The current `main.m` run creates no simulation figure because the initial
state plotting block is guarded by `if(0)`.

Plotting code is available in:

- `plot_dat_proc/threedplott.m` for vector-state figures or movies; and
- `plot_dat_proc/data_proc.m` for time traces and other analyses.

`threedplott.m` depends on an `arrow3` function that is not present in the
tracked repository. `main.m` also contains a disabled, machine-specific
`addpath` for that dependency.

## 5. Existing small benchmark cases

No current benchmark script, expected-output MAT-file, test, or documented
small benchmark case was found.

Git history shows that `main.m` previously used small lattices:

- `3 x 5` in commits including `bcad9bc` and earlier;
- `3 x 4` in later dipole/debug work; and
- `4 x 5` through commit `fe67a15`.

These are historical parameter settings, not reproducible baseline cases:
there are no accompanying fixed atom distributions, expected outputs,
tolerances, or benchmark instructions. Whether any of them was intended as a
general ASTRO benchmark is **UNCLEAR**.

## 6. Plotting convention for width and length

Use `plot_dat_proc/threedplott.m` as the plotting reference for preserving
the current physical convention:

- spin arrays are indexed as `(W, L, time)`;
- matrix rows / first dimension are width (`natomW`);
- matrix columns / second dimension are length (`natomL`);
- plotting maps length to the x-axis and width to the y-axis;
- x limits use `natomL`, and y limits use `natomW`; and
- the width index is reversed in the plotted position to preserve the
  established visual orientation.

The axis labels in `threedplott.m` are generic (`x axis`, `y axis`) rather
than explicit `length` and `width`, but its indexing and limits are the
clearest current implementation of the W/L convention. New benchmark plots
should follow that mapping rather than applying a plain matrix transpose.

Because `plot_dat_proc/threedplott.m` is currently untracked, whether it is
the intended maintained plotting API is **UNCLEAR**.

## 7. Most likely MATLAB command

From the repository root, the command most likely to run the current example
is:

```bash
matlab -batch "run('main.m')"
```

This command was not executed during the audit because it would run the full
GPU simulation and create `final.mat`, outside the requested read-only
benchmark investigation.

Expected prerequisites are MATLAB with Parallel Computing Toolbox support, a
compatible NVIDIA GPU, and a writable current directory.

## 8. Risks and unclear points

1. **No deterministic atom distribution.** `rng('shuffle')` makes every
   ordinary run different. A fixed `atomtype_` is safer than relying only on
   a MATLAB RNG seed across versions.
2. **No isolated parameter interface.** Parameters are embedded in
   `main.m`, which starts with `clear all`; a future benchmark needs a
   controlled entry point without changing the physics path.
3. **Full-workspace output is unstable.** `save('final.mat')` captures many
   implementation details and derived temporaries. A future baseline should
   define an explicit output schema.
4. **GPU numerical reproducibility.** Exact bitwise agreement across GPU
   models, CUDA stacks, MATLAB releases, or driver versions is **UNCLEAR**.
   Numerical tolerances and environment metadata will be required.
5. **Time vector alignment.** Each GPU chunk allocates `gpusteps` states but
   advances only `gpusteps-1` updates, while `runtime` and `t` are calculated
   from `gpusteps`. Chunk boundaries also save the carried state again.
   The exact physical timestamp corresponding to each saved spin plane is
   **UNCLEAR** and should be resolved before freezing expected trajectories.
6. **Composition naming.** The model distinguishes Gd from a combined FeCo
   TM species. The Co fraction is **UNCLEAR**.
7. **Magnetization contract.** Total and sublattice magnetizations are not
   produced by the simulation entry point. The post-processing formulas are
   not established as validated benchmark observables.
8. **No energy observable.** A physics-conservation or energy-based baseline
   cannot be made from current outputs without additional implementation.
9. **Plot dependency and status.** The available W/L-aware plotting function
   is untracked and requires missing `arrow3`; figure reproduction is not
   currently self-contained.
10. **Current run is not minimal.** The `20 x 30`, 10,000-nominal-step
    example is small enough for normal use but larger and slower than needed
    for a quick regression benchmark. The historical `3 x 5` case provides
    precedent, but its intended physics coverage is **UNCLEAR**.
11. **Historical cases are not references.** Small sizes in Git history do
    not include stored expected results or proof that they exercise all
    current fields correctly.
12. **Dipole modes.** Dipole mode 2 is explicitly marked as potentially
    wrong in source comments. The current baseline candidate avoids this by
    using `dipolemode=0`.
