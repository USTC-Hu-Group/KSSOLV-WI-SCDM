# WI tests in the KSSOLV branch

From this repository root:

```matlab
KSSOLV_startup;
summary = run_wi_tests;
```

The root launcher selects `tests/wi`; `run_release_tests` executes the original numerical assertions and a real PAO-library lookup against this tree's bundled UPF data. No separate setup or external solver-root setting is used. Scratch output is under `tests/wi/_work` or MATLAB temporary directories.

Coverage includes shifted/odd/even meshes and Fourier coefficients, band-error formulas and reference trimming, versioned metric caches/summaries, provenance success/failure records, 42 transform/inverse probes, three scalar/spinor MMN tests and eighteen fitting-policy checks. The fitting tests require Optimization Toolbox and are explicitly skipped if unavailable. The real-UPF lookup verifies that loaded PBE `ppvar.funct` overrides a PZ calculation label; it is count-only and does not certify full material projectability matrices.

The original test assertions/tolerances/data were preserved. Changes only adapt root discovery, `src/WI`, case-local scratch/fixture paths and the bundled startup. Test-only PAO/solver replacements live in `fixtures`; neither these nor all WI cases are on the normal startup path. Never recursively add this whole repository to the path.

This suite does not execute the entire original `test` collection, full material SCF/NSCF or Wannier90. See [validation evidence](../../docs/wi/VALIDATION_RESULTS.json) and the [release audit](../../RELEASE_AUDIT.md).

For source-package integrity, run `python tools/audit_branch.py` before publishing a clean source folder. That command does not clear license/PP redistribution requirements.
