# KSSOLV with WI-SCDM

This research branch integrates the WI-SCDM implementation into the supplied SOC-capable KSSOLV source tree. It includes KSSOLV solvers/classes, the supplied pseudopotential libraries, SCDM subspace construction, direct WI-SCDM/WI-SCDM-f interpolation, and an optional external Wannier90 interface.

**Local release candidate:** startup and numerical validation are recorded in [RELEASE_AUDIT.md](RELEASE_AUDIT.md). Public distribution still requires resolution of WI source-license provenance and review of bundled third-party/pseudopotential terms. No new project-wide license is assigned here. The original KSSOLV [README](README), source headers and bundled license files are preserved.

## Setup and quick start

Start MATLAB in this repository's root and use the ordinary KSSOLV startup:

```matlab
KSSOLV_startup;
report = run_wi_example('analytic');
```

There is no separate WI setup script, external KSSOLV checkout, or dependency-root environment variable to configure. `KSSOLV_startup` adds this tree's `src/WI` to the MATLAB path and runs its metric/path self-check. Case scripts and test fixtures are not added recursively. MATLAB path changes are session-local.

The analytic example generates small scalar and spinor states, constructs SCDM projections, runs both direct interpolators, and exercises overlap/Cube export. It checks orthogonality and input-mesh node energies. It does not perform DFT or execute Wannier90; off-mesh physical accuracy is not asserted. Outputs are under `example/wi/analytic/_output/`. See [example notes](example/wi/analytic/README.md) for expected complex-Cube and transformed-energy clipping warnings.

Existing KSSOLV APIs remain available after startup. For example, use `Atom`, `Crystal`, `setksopt`, `scf` and `nscf4c` as in the original KSSOLV examples. WI functions are called directly:

```matlab
% With compatible uniform-mesh KSSOLV states already computed:
% write_UNK_files(BX, crystal);
% write_input_files(crystal, Eigvals, [], 'erfc');
% [wan_r, Amn] = scdm_m('scdm.in', 'scdm.amn');
```

Run actual material protocols in dedicated case directories; the UNK writer replaces old case-local UNK files. Do not recursively add the whole repository with `genpath`: this would load test replacements and conflicting case stage names.

## Features

- Isolated/erfc/Gaussian scalar and spinor SCDM; multi-k projections, common QRCP selection and SVD orthogonalization.
- UPF PAO counting/projectabilities and automatic target/filter parameters, or manual filter/target choices.
- WI-SCDM and WI-SCDM-f on complete uniform meshes, including shifted/even meshes.
- Wannier90 input/overlap export, external MLWF execution and HR interpolation.
- Versioned band-error metrics, provenance records, orbital/spread and Hamiltonian-locality diagnostics.

WI projectabilities use an identity overlap metric and do not implement PAW/ultrasoft augmentation. The scalar/SOC protocols use the supplied norm-conserving ONCV files. Other KSSOLV capabilities and pseudopotential libraries are retained from the provided base, without new validation claims for every upstream module.

## Requirements

- MATLAB. Integration validation uses R2023a on Windows; no minimum release is established for this combined branch.
- Optimization Toolbox (`lsqcurvefit`) for automatic filter fitting; not required by the manual synthetic example or existing-file interpolation.
- The bundled KSSOLV tree/data supply the material API and PAO helpers. GPU/other optional upstream paths may need their original toolbox/platform dependencies.
- A separate Wannier90 executable for the MLWF stage only. Set `W90_COMMAND` for example drivers, or pass `W90Command` directly to `run_wannier90`. No executable is bundled.
- Python 3 standard library, optional, for `python tools/audit_branch.py`.

The supplied base had a platform-specific Linux `pbexsrF.mexa64` containing an author-specific build path. The release omits that binary and retains `pbexsrF.F90`; compile a local MEX before using the optional screened-exchange path. The documented WI demo/PBE preflights do not invoke that kernel. See [DEPENDENCIES.md](docs/wi/DEPENDENCIES.md).

## Material examples

Run each stage in a **fresh MATLAB process**, from the repository root:

```text
matlab -batch "run_wi_example('si','check')"
matlab -batch "run_wi_example('si','scf')"
matlab -batch "run_wi_example('si','nscf')"
matlab -batch "run_wi_example('si','scdm')"
matlab -batch "run_wi_example('si','wlmf')"
```

Replace `'si'` with `'soc'` for AuMn. `check` performs input/PP/PAO preflight only. The first four stages support direct interpolation; `wlmf` additionally requires Wannier90. No precomputed material checkpoints are bundled. The AuMn protocol has an inherited estimate of 46.85 GiB of UNK text and is not a small quick start.

See [Si](example/wi/si_scalar/README.md) and [AuMn/SOC](example/wi/soc_aumn/README.md). The SOC launcher enables the existing compatibility adapter before class construction/checkpoint loading. Manual SOC workflows should likewise call `enable_soc_blochham_compat` immediately after startup and before BlochHam use.

## Tests

```matlab
KSSOLV_startup;
summary = run_wi_tests;
```

This runs the existing WI numerical regressions and real bundled-UPF library lookup. Optimization-dependent checks are explicitly skipped when the toolbox is unavailable. Test fixtures stay isolated under `tests/wi/fixtures`. This suite does not run all original KSSOLV tests or full production DFT/MLWF calculations. See [test guide](tests/wi/README.md).

## Structure

```text
KSSOLV_startup.m       Ordinary KSSOLV startup with native WI path integration
run_wi_example.m      Select one bundled WI example
run_wi_tests.m        Run isolated WI numerical regressions
src/WI/               27 unchanged WI numerical source files
src/<other modules>/  Supplied KSSOLV classes, solvers and numerical modules
example/wi/           Analytic, Si and AuMn/SOC protocols
example/, test/       Original KSSOLV examples/tests retained
tests/wi/             WI tests and isolated replacement fixtures
ppdata/, external/    Supplied PP libraries and third-party components
util/, structure/     Original KSSOLV utilities and structures
doc/                  Original KSSOLV documentation
docs/wi/              WI guides, copied-source manifest and validation evidence
```

[CODE_STRUCTURE.md](docs/wi/CODE_STRUCTURE.md) describes the numerical chain. [SOURCE_MANIFEST.json](docs/wi/SOURCE_MANIFEST.json) fingerprints each included source against its supplied component. This identifies the packaged snapshot even though its upstream Git commit has not been established. No Git history or official upstream branch relationship is fabricated.

## Citation

Final WI publication/software identifiers are pending. Replace this paragraph with the approved author/title/journal/year/DOI and software identifier when assigned. Cite the relevant KSSOLV, SCDM, transformed-interpolation and Wannier90 publications separately as appropriate. No WI DOI or remote repository URL is invented here.

## License

Read [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). KSSOLV files carry upstream notices; bundled components include their own licenses, including GPLv3 HANSO/NLCG and the export_fig BSD-style notice. The WI translation-origin GPL statement conflicts with a manuscript BSD field and has not been resolved by merging folders. Pseudopotential redistribution terms also require review. No blanket MIT/BSD license or claim of sole WI authorship is applied to this branch.
