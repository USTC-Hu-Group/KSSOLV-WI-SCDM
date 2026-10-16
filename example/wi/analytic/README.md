# Analytic scalar/spinor demo

From the repository root:

```matlab
KSSOLV_startup;
report = run_wi_example('analytic');
```

No external KSSOLV checkout or nonlinear-fitting toolbox is needed. The bundled solver startup sets all paths. The example uses seed 17, four real-space grid points, three k points, three bands, a two-orbital target and manual erfc parameters `mu=0`, `sigma=0.5` eV. It exercises the unchanged SCDM, both direct interpolators, MMN/WIN writers, Cubes and spread estimates. Its local helper writes synthetic text UNK directly; it does not test the material inverse-FFT exporter.

Outputs are created under `_output/<unique-run>/scalar` and `spinor`, including `demo_results.mat` and `summary.json`. Assertions require orthogonality below `1e-10` and input-node energies within `1e-8` eV. Current combined-branch validation passed with MATLAB R2023a; exact records are in [VALIDATION_RESULTS.json](../../../docs/wi/VALIDATION_RESULTS.json).

This is not a material benchmark. The one-neighbour `toy.nnkp` is an explicit overlap test fixture, not a complete Wannier90 finite-difference set; no external Wannier90 execution occurs. Scalar Cubes retain only phase-aligned real parts and can warn about residual imaginary weight. The chosen synthetic spinor f-interpolation can produce positive transformed off-mesh eigenvalues, clipped to epsilon by the inherited inverse (about 0.115 eV maximum in the run). Off-mesh physical accuracy and Cube spread fidelity are not asserted; complex `wan_r` is preserved in MATLAB output.
