# Changes relative to supplied component copies

All 27 numerical files under `src/WI` are byte-identical to the WI source release. Solver algorithms, class constructors, pseudopotential contents, numerical settings and existing test assertions were not rewritten. Root-discovery/path/diagnostic/document changes are listed below. UTF-8/LF normalization accompanies edited text copies.

| Component / destination | Changed release line spans | Reason / effect |
| --- | --- | --- |
| `kssolv` / `.gitignore` | 2-48 | Merge output/privacy ignores and retain the intentional upstream GUI FIG asset. |
| `kssolv` / `doc/htmldoc/source/Start/scf_setup.rst` | 18-19 | Replace illustrative upstream installation/user paths with portable placeholders; upstream notices retained. Historical prebuilt manual PDF is retained without recompilation. |
| `kssolv` / `doc/texdoc/kssolvug.tex` | 200, 204 | Replace illustrative upstream installation/user paths with portable placeholders; upstream notices retained. Historical prebuilt manual PDF is retained without recompilation. |
| `kssolv` / `external/export_fig/README.md` | 156 | Replace illustrative upstream installation/user paths with portable placeholders; upstream notices retained. Historical prebuilt manual PDF is retained without recompilation. |
| `kssolv` / `KSSOLV_startup.m` | 31-35 | Append native WI path precedence and metric/path self-check; preserve original CPU startup and upstream header. |
| `kssolv` / `util/outpwscf.m` | 35 | Replace author-specific QE pseudo-directory output default with case-relative ./pseudo/; no WI numerical effect. |
| `wi_release` / `docs/wi/CODE_STRUCTURE.md` | 3, 9, 86 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `docs/wi/ERROR_METRICS_VITALE2020.md` | 42 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `docs/wi/MODULE_SELECTION.csv` | 1-36 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `example/wi/analytic/README.md` | 1, 6-7, 10, 12, 14 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `example/wi/analytic/run_analytic_demo.m` | 6, 9 | Use bundled startup and new directory depth, or replace obsolete root-variable diagnostic hints. Computation settings unchanged. |
| `wi_release` / `example/wi/si_scalar/README.md` | 1, 3, 5, 7-12, 15, 17, 19, 21 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `example/wi/si_scalar/run_si_pipeline.m` | 5, 7-8 | Use bundled startup and new directory depth, or replace obsolete root-variable diagnostic hints. Computation settings unchanged. |
| `wi_release` / `example/wi/soc_aumn/README.md` | 1, 3, 5, 7, 10-14, 17, 19, 21 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `example/wi/soc_aumn/run_scdm.m` | 78 | Use bundled startup and new directory depth, or replace obsolete root-variable diagnostic hints. Computation settings unchanged. |
| `wi_release` / `example/wi/soc_aumn/run_soc_pipeline.m` | 5, 7 | Use bundled startup and new directory depth, or replace obsolete root-variable diagnostic hints. Computation settings unchanged. |
| `wi_release` / `example/wi/soc_aumn/run_wlmf.m` | 65 | Use bundled startup and new directory depth, or replace obsolete root-variable diagnostic hints. Computation settings unchanged. |
| `wi_release` / `tests/wi/audit_wi_scdm_f_inverse.m` | 3, 6-7 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/README.md` | 1, 3, 6-7, 10, 12, 14, 16, 18 | Rewrite/adapt branch documentation, module destinations, commands or relative links to the integrated layout. |
| `wi_release` / `tests/wi/run_release_tests.m` | 1-3, 6-7, 9, 19, 28-29 | Use bundled startup, nested test paths and unconditional real bundled-UPF lookup; existing numerical tests unchanged. |
| `wi_release` / `tests/wi/test_calculate_band_error_modes.m` | 3-4 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_calculate_err_trims_reference_bands.m` | 4-5 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_compute_pao_projectability_ppvar_funct.m` | 4-5, 7, 8 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_scdm_fit_policy.m` | 3, 6-7, 11, 67, 89 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_scdm_provenance.m` | 3-4 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_scdm_shifted_kgrid.m` | 4-7 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_vitale_error_pipeline.m` | 3-4 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_vitale_error_summary.m` | 3-4 | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |
| `wi_release` / `tests/wi/test_write_w90_mmn_from_unk.m` | Newline/encoding only | Adapt repository/src/WI/scratch/fixture paths and bundled startup. Assertions and numerical values unchanged; any remaining differences are text normalization. |

New branch-owned packaging files: root README, `run_wi_example.m`, `run_wi_tests.m`, dependency/test/example guides, third-party notices, branch integrity tool, source/validation/change manifests and release audit. No standalone WI setup script is included.

Excluded inherited files (not removed from the original component):

| File | Reason |
| --- | --- |
| `src/TDDFT/out.mat` | Unreferenced generated TDDFT result |
| `src/SCF_DCM/Energy/getEewald.asv` | MATLAB autosave |
| `src/SCF_DCM/Spin/iterinit/.genX0_random.m.swp` | Editor swap |
| `src/BasisClass/@BlochHam/BlochHam.m.bak_soc_nonmag_20260910_160936_405055` | Historical class backup |
| `src/SCF_DCM/Potential/Vxc/VGCExchange/pbexsrF.mexa64` | Optional platform-specific MEX contains author-specific build path; Fortran source retained, rebuild locally for screened exchange. |
