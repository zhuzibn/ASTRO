# ASTRO
ASTRO (atom spintronic dynamics) is an atomistic model developed by COSMAL group for the simulation of spin dynamics in 2D ferrimagnet GdFeCo. This model can be modified to support simulating ferromagnet, antiferromagnet in 2D or 3D structure.
If you used this environment for your experiments or found it helpful, consider citing the following papers:
1D FiM: [Phys. Rev. Appl. 13, 034040 (2020)]
2D FiM: [Phys. Rev. B 106, 184419 (2022)]
3D FiM: [Appl. Phys. Lett. 124, 012405 (2024)]
AFM: [Phys. Rev. B 109, 134433 (2024)]

## Verification

Production simulation defaults are defined in `astro_default_config.m` and
consumed by both `main.m` and the benchmark runner. Benchmark runs override
only the fixed atom distribution and output path; quick mode also overrides
the lattice size.

Deterministic atom-distribution validation and initial spin construction are
implemented in `astro_validate_atom_distribution.m` and
`astro_initial_spin_state.m`. `systemgeneration.m` still owns random,
load/save, and restart policy, while the pure initializer preserves the
current RE/Gd and TM/FeCo sign conventions exactly.

Run the deterministic invariant smoke checks from the repository root:

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
```

The smoke command runs the quick benchmark into a temporary candidate
directory and exits nonzero if the fixed atom distribution, trajectory shape,
finite spin values, spin normalization, or saved-time contract regresses.
It also checks that unsupported solver selections `rk4=0` and `rk4=2` fail
with the documented `ASTRO:UnsupportedSolver` error before integration, and
that deterministic initialization matches the previous fixed-input formula
exactly.

Run deterministic boundary and field-term checks from the repository root:

```bash
matlab -batch "run('tests/run_field_tests.m')"
```

The field command checks `bc=1` non-periodic and `bc=0` periodic neighbor
behavior on small fixed atom matrices. It also verifies exchange,
anisotropy, DMI, and external-field terms separately against hand-derived
expected values while leaving thermal and dipole terms out of scope.
