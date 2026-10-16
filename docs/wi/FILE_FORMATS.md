# File, coordinate and energy contracts

## Arrays and units

| Object | Shape / convention |
| --- | --- |
| KSSOLV `Eigvals` | `[nbnd,nk]`, Hartree; material scripts reshape solver output explicitly. |
| `crystal.efermi` | Hartree. Interpolators convert with `27.211386246` eV/Ha and subtract this same DFT reference. |
| Finalized `crystal.kpts` | Rows of Cartesian reciprocal coordinates, Bohr inverse. Reciprocal lattice rows are `2*pi*inv(C).'` for cell-row matrix C in Bohr. |
| SCDM/input and query k points | Reciprocal fractional coordinates, three columns. Complete uniform construction mesh; arbitrary query paths. |
| `Amn{k}` | `[nbnd,num_wann]`, complex. One per input k point. |
| `wan_r` | `[2*Ngrid,num_wann]`, complex, spin up then spin down. Scalar data have a zero lower block. |
| Direct/HR interpolated bands | `[nk_path,num_wann]`, eV relative to DFT Fermi energy. |
| `calculate_band_error` inputs | `[band,kpoint]`, eV relative to the **same** DFT Fermi energy, matched bands and k points. |
| Direct `H_R` | `[N1,N2,N3,num_wann,num_wann]`, centered integer-R convention and mesh-origin phase. WI-SCDM-f stores transformed coefficients. |

Manual filter `mu/sigma` and transform `epsilon/a` are absolute eV, matching `scdm.eig`; metric `Nu/Tau` are eV relative to DFT Fermi energy. Confusing these references changes the calculation.

## Files produced and consumed

| File | Producer -> consumer | Contract |
| --- | --- | --- |
| `tmp/UNK%05d.1` | `write_UNK_files` -> SCDM/PAO/MMN | Formatted text, header `nx ny nz ik nbnd`; cell-periodic scalar `u_nk`, pairs `Re Im`. |
| `tmp/UNK%05d.NC` | Same | Same header; interleaved `Re_up Im_up Re_down Im_down`, grid rows grouped by band; full-spinor normalization. |
| `scdm.in` | `write_input_files` -> SCDM/interpolators/export | First line `nbnd nk num_wann`; second line eigenvalue filename, filter and applicable parameters; next nk lines fractional k points. Trailing metadata includes `KGRID`, `CELL_BOHR` and atoms. Use the writer rather than hand-editing. |
| `scdm.eig` | `write_input_files` / `write_w90_eig` -> SCDM/W90 | One-based band and k indices followed by absolute eV; band index runs fastest. |
| `scdm.amn` | `scdm_m` -> direct interpolation/W90 | Comment, `nbnd nk num_wann`, then `band wannier k Re Im`, one-based. Active AMN writer is internal to `scdm_m`. |
| `seed.win` | `write_w90_win` -> W90 | Cell, atom/k coordinates, mesh, counts and optimization settings. Material wrappers use finalized KSSOLV objects. |
| `seed.nnkp` | W90 `-pp` -> `write_w90_mmn` | Neighbour index and reciprocal shift G; must match mesh, bands and UNK states. |
| `seed.mmn` | `write_w90_mmn` -> W90 | Neighbour overlaps, reciprocal-shift phases and both spinor components. |
| `seed_hr.dat`, optional WS-vector data | W90 -> `compare_w90_bands` | Physical Hamiltonian and degeneracy/WS information; local helpers perform HR interpolation. |
| `scdm_*.cube` | `write_wannier_cubes` -> spread diagnostic | Scalar phase-aligned real part or spinor total amplitude; not lossless complex orbitals. |
| `scdm_provenance/` | `save_scdm_provenance` -> audit | JSON/MAT plus source and PP snapshots; runtime-only, may include upstream code and local paths. |
| `err.txt` | `calculate_err` -> metric cache/summary | Append-only records: method, formula version, mode, `Nu/Tau`, band selection and dimensions. |

UNK here means formatted text handled by this project. Do not substitute binary UNK without conversion. `scdm_m` checks the input directory and UNK locations with local readers; PAO construction expects case-local `tmp`. Material launchers change into their own case before writing files.

## Material checkpoint chain

`scf_checkpoint.mat` supplies density and converged state to `run_nscf`. `nscf_checkpoint.mat` supplies the uniform-mesh Crystal, Basis and eigenvalues; `nscf_path_checkpoint.mat` supplies independent DFT path energies. SCDM writes projection/orbital and direct-interpolation checkpoints; MLWF writes a separate checkpoint and Wannier90 output directory. These class checkpoints depend on the KSSOLV implementation and are not bundled.
