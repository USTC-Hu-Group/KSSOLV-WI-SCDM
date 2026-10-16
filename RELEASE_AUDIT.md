# KSSOLV/WI combined-branch release audit

Date: 2026-10-15. Local integration only; no Git initialization, branch history, remote repository, commit or push was created.

## Result and layout

The supplied KSSOLV snapshot and WI source release are integrated into one KSSOLV-style source tree. Ordinary `KSSOLV_startup` exposes the WI module. The standalone WI setup file is absent, and published MATLAB launchers/tests do not read external KSSOLV/WI-root environment variables.

```text
KSSOLV_WI_SCDM_Release/
|-- KSSOLV_startup.m          # native startup + WI path/metric checks
|-- run_wi_example.m          # analytic / si / soc, one case at a time
|-- run_wi_tests.m
|-- README.md                # current branch guide
|-- README                   # preserved original KSSOLV guide
|-- THIRD_PARTY_NOTICES.md
|-- src/
|   |-- WI/                  # 27 unchanged WI numerical functions
|   `-- <KSSOLV modules>/    # original solver/class layout retained
|-- example/wi/              # analytic, scalar Si, AuMn/SOC
|-- example/, test/          # original KSSOLV examples/tests retained
|-- tests/wi/                # regressions and isolated fixtures
|-- util/, structure/
|-- ppdata/                  # supplied libraries, original data unchanged
|-- external/                # existing source/assets/license files
|-- doc/                     # original upstream manuals/sources
|-- docs/wi/                 # technical, provenance and validation records
`-- tools/audit_branch.py
```

The KSSOLV module tree was preserved instead of reducing a software branch to WI-only files. WI's previous dependency-based 27-function selection remains, including PAOs, both direct methods, scalar/spinor export, Wannier90 interfaces, current error metrics and numerical spread/locality diagnostics. Test replacements are in `tests/wi/fixtures`, outside the `test` directory recursively added by the upstream startup. Case directories are selected individually by the root launcher to avoid duplicate `run_scf`/`run_nscf` names.

## Original protection and provenance

The two supplied component directories were only read. Before copying, all 679 KSSOLV files and all 74 WI-release files were recorded with SHA-256, size and modification time. Final verification covers all 753 originals. The combined snapshot is identified by [SOURCE_MANIFEST.json](docs/wi/SOURCE_MANIFEST.json), containing 739 included copied-file entries and source/release hashes. It does not fabricate an upstream commit or official fork relationship.

The key BlochHam, NSCF, PpData, library selector and Si/Au/Mn data match the previously validated local SOC KSSOLV copy. Each of the 27 WI numerical sources remains byte-identical to the earlier WI package. Startup/path integration changes the locations and selected runtime tree, not the mathematical implementation.

## Copy changes and exclusions

All changed copied files, line spans, reasons and impacts are listed in [SOURCE_CHANGES.md](docs/wi/SOURCE_CHANGES.md). Changes comprise native startup integration, portable example/test root discovery, removal of obsolete environment-root dependencies/hints, documentation/relative-link updates, output ignores and portable legacy QE-export/document examples. Original source copyright/license headers remain intact.

Five inherited artifacts are excluded from this copy only: a MATLAB autosave, an editor swap, a historical BlochHam backup, an unreferenced TDDFT `out.mat`, and the optional Linux `pbexsrF.mexa64`. The MEX contained an author-specific build path; its GPL-header Fortran source is retained for a local rebuild. The documented PBE WI checks do not select this optional screened-exchange kernel. It must be rebuilt/validated before using that upstream path. The upstream GUI FIG and Java class/source assets are retained because they are source/runtime assets rather than generated WI results; `.gitignore` explicitly retains the GUI FIG.

No original file was deleted, and no pseudopotential numerical data, solver/class algorithm, WI numerical file, physical setting or existing test assertion was rewritten. Edited upstream document source examples use portable placeholders; the prebuilt upstream PDF remains a historical manual and was not recompiled. Current installation commands are in the root README.

## Validation performed on the combined tree

Detailed evidence is in [VALIDATION_RESULTS.json](docs/wi/VALIDATION_RESULTS.json). The older [standalone record](docs/wi/STANDALONE_VALIDATION_RESULTS.json) is clearly labeled historical and is not substituted for current integration tests.

- MATLAB R2023a Windows, fresh default MATLAB path, external-root settings cleared: normal startup and resolved core/solver paths point inside this combined tree. No standalone setup function was found and fixtures were not on the normal startup path.
- Root analytic launcher passed scalar/spinor SCDM, both direct interpolators and MMN/WIN/Cube/spread operations. Orthogonality and input-node energies meet the existing assertions (about 1e-15 and 1e-13 eV respectively). Expected complex-scalar Cube and off-mesh transformed-energy clipping warnings remain documented; no off-mesh material-accuracy claim is made.
- Root test launcher passed **10 groups, zero skips**, including three MMN function tests, eighteen fit-policy checks with real Optimization Toolbox fitting, twelve shifted-grid/method combinations, 42 inverse probes, metric/provenance/summary tests and real bundled-UPF functional lookup.
- Si preflight passed against bundled classes/data: 216 SCF k points and eight PAOs. AuMn/SOC preflight passed in a separate fresh process: 64 SCF points and 72 PAOs with fully relativistic headers/valences and compatibility adapter. These are input checks, not SCF/NSCF runs.
- `checkcode` inspected 58 WI/launcher/modified-startup/utility MATLAB files: no parse errors, 17 nonfatal messages. Ten originate in the prior WI/test files, five in the legacy QE-export utility, one in the original startup unused assignment and one in the new `feval` dispatcher. Locations/IDs are recorded. No algorithm/style refactoring was applied.
- Static required-file analysis was refreshed after the optional binary exclusion: 159 resolved custom-code files, all inside the combined tree, and no custom code outside it/MATLAB. Dynamic PP, provenance, SOC class and external-command lookups were reviewed separately. Optional screened-exchange MEX rebuilding remains explicit rather than being treated as statically certified.
- Source-package checker verifies copied-file hashes, unchanged WI core, relative Markdown links, obvious private paths/tokens, runtime-data exclusions, optional binary provenance and absence of standalone-setup/external-root calls. Final result is recorded below.

Raw logs, temporary class/source/PP snapshots, synthetic results and test outputs were removed after recording sanitized evidence. All test processes used new-tree case/scratch directories; original component directories were not used as output locations.

## Validation and publication limits

Full material SCF/NSCF, real-material projectability matrices/fitting, external Wannier90 execution, the complete inherited KSSOLV test collection and screened-exchange/cross-platform rebuilds were not executed here. The optional GPU dispatch still has its upstream toolbox requirements. No Wannier90 executable is bundled or newly verified. Do not equate source/path/overlap tests with production physical validation.

The combined tree removes the need to obtain a separate matching local KSSOLV installation for the documented checks. Public-upstream version identification is still a provenance task; the bundled hashes identify the exact provided snapshot meanwhile.

**Public distribution is not yet cleared.** The original WI README identifies translated GPL source while a manuscript field says BSD. The bundled tree also includes GPLv3 HANSO/NLCG, a GPL-header screened-exchange source, export_fig's own notice and PP data whose redistribution provenance requires review. All available notices/texts were retained, but no blanket project license was invented. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Merging folders does not resolve these rights questions. Final citation/software identifiers and production validation should be completed before paper-reproducibility claims.

The package/privacy checker is bounded, not a guarantee about every historical upstream document or legal obligation. The unchanged WI path validator retains a legacy error hint mentioning `WI_ROOT`; it does not read that setting. In this branch the hint means checking the active `src/WI` path.

## Final structural check

`python tools/audit_branch.py`: **PASS**. Checked 749 delivered source/input/document/asset files, 739 copied-file hashes and 27 unchanged WI numerical files. No detected credential files/tokens, author-specific runtime paths, broken relative Markdown links, standalone setup calls, external-root environment reads or generated/scratch files remain. All resolved static dependencies recorded for the documented closure exist inside this tree. The package is approximately 51.3 MiB; no individual source/input/asset exceeds 1 MiB. Available third-party notices and source/data provenance still require the separate public-distribution review stated above. Original-preservation verification: all 753 source-component files retained their initial SHA-256, size and modification time.
