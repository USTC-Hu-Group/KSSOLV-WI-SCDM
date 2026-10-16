# Scalar diamond Si

The existing protocol uses a two-atom primitive diamond cell with lattice constant 10.263 Bohr, PBE, cutoff 30 Ha, 300 K Gaussian smearing and 12 bands. SCF uses 6x6x6 k points, uniform NSCF uses 8x8x8 (512 points), and the G-X-W-K-G-L comparison path has 151 points. The bundled scalar Si UPF gives eight PAOs per cell. Automatic erfc fitting sets the construction parameters; both direct methods use the same SCDM basis.

From the repository root, run one command per fresh MATLAB process:

```text
matlab -batch "run_wi_example('si','check')"
matlab -batch "run_wi_example('si','scf')"
matlab -batch "run_wi_example('si','nscf')"
matlab -batch "run_wi_example('si','scdm')"
matlab -batch "run_wi_example('si','wlmf')"
```

KSSOLV startup and case paths are selected automatically by the launcher. No separate WI setup or external KSSOLV root setting is needed. Automatic fitting needs Optimization Toolbox. The optional last stage requires an external Wannier90 command, supplied by `W90_COMMAND` or the default `wannier90.x`; see [dependencies](../../../docs/wi/DEPENDENCIES.md).

`check` validates input/PAO count without SCF. SCF must converge before NSCF. NSCF writes both uniform states and independent DFT path energies; `scdm` writes projections, both direct band arrays, errors and spread estimates. `wlmf` runs external gauge refinement and HR comparison against that same reference. Stop at `scdm` for the direct methods.

All outputs stay under this case: SCF/NSCF/path/SCDM/comparison/MLWF checkpoints, `tmp/UNK*`, SCDM input/eigenvalue/projection/Cube files, `err.txt`, `logs`, `results` and `wannier90_post`. Results include DFT/per-method band CSVs, `si_error_summary.csv` and spread estimates. Runtime products are not delivered.

The **combined branch** passed preflight using its bundled PP and classes: eight PAOs and 216 SCF points. Full SCF/512-point NSCF/material fitting and Wannier90 were not rerun in this integration. The NSCF driver explicitly checks 512 points; changing only the mesh is not a supported small-case reduction. Use the analytic demo for a minimal run.
