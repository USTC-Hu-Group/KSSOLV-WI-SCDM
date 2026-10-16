# Bundled-branch dependencies

## Included snapshot

The supplied `kssolv-wi-scdm` tree is bundled with its `src`, `util`, `structure`, `example`, `test`, `doc`, `external` and `ppdata` directories. WI numerical functions are integrated under `src/WI`. A separate KSSOLV installation is neither required nor selected by environment variables. Use `KSSOLV_startup` from this root.

The key BlochHam/NSCF/PpData/library-selector files and the Si/Au/Mn UPF inputs were compared bytewise with the previously used local SOC KSSOLV copy and match it. The exact public upstream branch/commit is still unconfirmed. `SOURCE_MANIFEST.json` records the supplied snapshot hashes; it does not claim an official upstream revision.

`KSSOLV_startup` selects CPU through the existing `switchver` routine, adds the original module paths, then explicitly prefers `src/WI` and validates metric paths. WI cases and fixture replacements are excluded from recursive startup. Do not add the whole repository recursively.

## Runtime requirements

| Dependency | Scope |
| --- | --- |
| MATLAB | All numerical code; R2023a Windows is the integration-check environment, not an established minimum version. |
| Optimization Toolbox | `lsqcurvefit` for automatic filter fitting; manual target/filter interpolation and analytic demo avoid it. |
| Parallel Computing Toolbox/GPU | Optional upstream GPU eigensolver alternatives. Not required by the CPU synthetic demonstration or direct-file interpolation. |
| Wannier90 executable | External preprocessing/MLWF only. Not included; actual executable execution is not certified by overlap tests. |
| Fortran/MEX toolchain | Rebuilding platform-specific optional upstream kernels. The inherited Linux binary was excluded because it embeds an author-specific build path. `pbexsrF.F90` is retained; build a local MEX for the optional screened-exchange path. This build was not executed during integration. |
| Export/rendering tools | Optional legacy export_fig functionality may require Ghostscript/Xpdf/Java under its original requirements. Generic Program Files discovery defaults are upstream software-search paths, not author-specific settings. |
| Python 3 | Optional standard-library package audit only. No Python transform is required by WI. |

The original upstream documentation in `doc` and source headers remain relevant to KSSOLV modules not exercised by WI tests. Presence of their source is not a claim that every module/build path works on every platform.

## Pseudopotentials

WI examples use these bundled files, relative to this repository root:

| File | Intended use |
| --- | --- |
| `ppdata/default/Si_ONCV_PBE-1.3.upf` | Scalar diamond Si; four PAOs per atom and eight per cell. |
| `ppdata/default_FR/Au_ONCV_PBE_FR-1.3.upf` | Fully relativistic Au, valence 19 and 18 spatial PAOs. |
| `ppdata/default_FR/Mn_ONCV_PBE_FR-1.3.upf` | Fully relativistic Mn, valence 15 and 18 spatial PAOs. |

SOC doubles the spatial PAO count, producing 72 for AuMn. Fingerprints are in `SOURCE_MANIFEST.json` and the case JSON. The PAO code locates `ppdata` from this tree's `kssolvpptype` path and prefers the functional actually loaded in `ppvar`. It reads `PP_PSWFC/PP_CHI.*` and uses `Ggrid`, `ylmr`, `bessel`, `simpson` from the bundled solver.

Projectabilities use an identity plane-wave metric without PAW/ultrasoft augmentation. The base also contains other PP libraries for upstream functionality; their inclusion does not extend WI's overlap treatment. Review PP source/distribution rights before public release; a checksum is not a license.

## SOC

Use `run_wi_example('soc',stage)` in a fresh process per stage. The existing `enable_soc_blochham_compat` adapter runs before object construction/loading, recognizes the known class layout and creates a temporary adapted class while preserving the upstream header. It does not edit this repository's solver class. For a manual SOC workflow, call it immediately after startup and before using BlochHam.

## Wannier90

Obtain an executable from the [official distribution](https://wannier.org/download/). Historical scripts referred to 3.1.0; no particular executable version is verified for this branch. Example launchers use `W90_COMMAND`; the lower-level API accepts `W90Command` or falls back to `WANNIER90`, then `wannier90.x`.

```matlab
setenv('W90_COMMAND','wannier90.x');
run_wi_example('si','wlmf'); % requires the completed preceding material stages
```

The interface runs `wannier90.x -pp seedname`, builds MMN from the resulting NN packet and UNK data, then runs `wannier90.x seedname`. Default `dis_num_iter=0` optimizes the SCDM subspace's gauge; changing it changes the scientific comparison. Windows paths with spaces require suitable shell quoting of the executable command.

## Validation boundaries

Current branch validation is recorded in `VALIDATION_RESULTS.json` and the root [audit](../../RELEASE_AUDIT.md). `STANDALONE_VALIDATION_RESULTS.json` is explicitly a **historical** record from the earlier WI-only package, not evidence that the merged branch was already tested then. Complete material SCF/NSCF, full material PAO fitting and external MLWF execution need separate production validation.

## Optional screened-exchange kernel build

Configure a MATLAB-supported Fortran compiler with `mex -setup FORTRAN`, then from `src/SCF_DCM/Potential/Vxc/VGCExchange` build `pbexsrF.F90` with the compiler/preprocessor settings appropriate to that MATLAB/platform. The source includes MATLAB `fintrf.h`; for a configured toolchain, `mex pbexsrF.F90` is a starting command, not a tested cross-platform build recipe. Keep the resulting platform binary local (ignored by Git). Full hybrid/screened-exchange calculations require validation after rebuilding. The PBE WI examples do not select this optional kernel.

The unchanged WI path validator contains a legacy diagnostic mentioning `WI_ROOT`. In this branch it means checking the active `src/WI` function path; the validator does not read that environment variable. Normal startup and launchers require neither external-root variable. Historical upstream manuals under `doc` are retained for reference; use this branch README for current installation and WI commands.
