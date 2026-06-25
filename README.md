# ASTRO
ASTRO (atom spintronic dynamics) is an atomistic model developed by COSMAL group for the simulation of spin dynamics in 2D ferrimagnet GdFeCo. This model can be modified to support simulating ferromagnet, antiferromagnet in 2D or 3D structure.
If you used this environment for your experiments or found it helpful, consider citing the following papers:
1D FiM: [Phys. Rev. Appl. 13, 034040 (2020)]
2D FiM: [Phys. Rev. B 106, 184419 (2022)]
3D FiM: [Appl. Phys. Lett. 124, 012405 (2024)]
AFM: [Phys. Rev. B 109, 134433 (2024)]

## Verification

Run the deterministic invariant smoke checks from the repository root:

```bash
matlab -batch "run('tests/run_smoke_tests.m')"
```

The smoke command runs the quick benchmark into a temporary candidate
directory and exits nonzero if the fixed atom distribution, trajectory shape,
finite spin values, spin normalization, or saved-time contract regresses.
