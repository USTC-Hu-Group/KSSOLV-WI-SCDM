# AuMn noncollinear SOC

The inherited two-atom AuMn input uses PBE, 40 Ha, `nspin=4`, noncollinearity/SOC, 34 valence electrons and initial Mn moment 3 along z. SCF is 4x4x4/42 bands; NSCF is 6x6x6/90 bands. The generic G-K1-K2-K3-G-K4 path has 121 points; these are not automatically assigned symmetry labels. Fully relativistic Au/Mn data give 72 spin-component PAOs.

`structure.xsf` is the reference structure in Angstrom; the stage script specifies the equivalent Bohr cell/coordinates. `case_metadata.json` records structure/PP hashes and inherited storage estimates; it is provenance rather than a runtime settings loader.

From the repository root, use a fresh MATLAB process per stage:

```text
matlab -batch "run_wi_example('soc','check')"
matlab -batch "run_wi_example('soc','scf')"
matlab -batch "run_wi_example('soc','nscf')"
matlab -batch "run_wi_example('soc','scdm')"
matlab -batch "run_wi_example('soc','wlmf')"
```

The launcher initializes this bundled KSSOLV, selects the case, and activates `enable_soc_blochham_compat` before class construction/loading. No separate WI setup or external KSSOLV-root configuration is needed. PP inputs come from the bundled `ppdata/default_FR`; fitting requires Optimization Toolbox. `wlmf` additionally requires external Wannier90 (`W90_COMMAND`) and requests orbital plotting. See [dependencies](../../../docs/wi/DEPENDENCIES.md).

`check` validates SOC APIs, fully relativistic PP headers/valences and PAO count without SCF. Subsequent stages need preceding converged checkpoints. Outputs include SCF/uniform/path/SCDM/comparison checkpoints, `tmp/UNK*.NC`, SCDM/overlap/Cube files, errors, logs and `wannier90_post`.

The combined branch's preflight passed with its own classes/PP data: 64 SCF points and 72 PAOs. Full SOC SCF/NSCF, magnetic convergence, material fitting and external MLWF were not repeated. The inherited estimate of 46.85 GiB text UNK is not a newly measured storage or memory requirement. This protocol is not the minimal demo and is not a newly certified magnetic ground state.
