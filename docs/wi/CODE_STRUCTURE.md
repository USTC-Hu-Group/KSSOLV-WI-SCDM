# Code structure and numerical workflow

This guide summarizes the inspected implementation, including internal helpers and runtime lookups. It supersedes prototype interfaces described in the original project README. The release keeps 27 existing MATLAB functions in the `src/WI` directory of the bundled KSSOLV branch; none of those files was numerically rewritten. Source fingerprints and the selection of all 35 original `src` files are recorded in `SOURCE_MANIFEST.json` and `MODULE_SELECTION.csv`.

## Entry points and modules

| Area | Functions | Responsibility |
| --- | --- | --- |
| Setup and validation | `KSSOLV_startup`, `validate_vitale_error_runtime`, `validate_soc_kssolv_runtime`, `enable_soc_blochham_compat` | Initialize bundled KSSOLV and prefer this WI source tree; check metric paths/formulas and SOC APIs. The SOC adapter creates a temporary class shadow before class loading. |
| Bloch-state export | `write_UNK_files` | Inverse FFT of KSSOLV Basis coefficients, normalized scalar/spinor cell-periodic states in text UNK files. |
| Inputs and parameter selection | `write_input_files`, `compute_pao_projectability`, `save_scdm_provenance` | Hartree/eV conversion, grid/structure metadata, PAO target count/projectability, optional filter fit and reproducibility records. |
| SCDM subspace | `scdm_m`, `read_scdm_uniform_kgrid` | Anchor selection, weighted QR with column pivoting, all-k projections, polar factor by SVD, AMN and real-space orbitals. UNK and AMN helpers are local functions inside `scdm_m`. |
| Direct interpolation | `compute_scdm_only_bands`, `compare_scdm_bands` | Projected Hamiltonian and grid-aware Fourier interpolation; the comparison wrapper also builds a path and plots bands. |
| Transformed interpolation | `wi_scdm_f_params`, `wi_scdm_f_transform`, `wi_scdm_f_inverse` | Energy transform selection/evaluation and step-limited Newton inverse, used inside direct interpolation. |
| Wannier90 interface | `run_wannier90`, `write_w90_win`, `write_w90_eig`, `write_w90_mmn`, `compare_w90_bands` | Prepare files, run preprocessing, generate overlaps, run MLWF optimization and read/interpolate HR output. HR readers/interpolators are local helpers in the comparison wrapper. |
| Error evaluation | `calculate_band_error`, `calculate_err`, `has_vitale_band_error`, `run_eta_only`, `collect_eta_summary` | Independent Fermi-Dirac error weights, common Fermi reference, versioned records, recalculation and matched summaries. |
| Real-space diagnostics | `write_wannier_cubes`, `calculate_scdm_spread`, `plot_scdm_locality` | Cube representation, orbital-center/spread estimates, numerical decay of Hamiltonian blocks. Locality is a computational diagnostic as well as a plotting routine. |

`write_w90_eig` is a separately callable export utility; `write_input_files` already writes `scdm.eig`. `run_wannier90` copies/renames that existing file. The older standalone `write_amn` and `readUNK` are not dependencies: the active routines contain their own scalar/spinor-aware serializers/readers.

## Main flow

```text
External KSSOLV: Atom / Crystal / setksopt -> scf -> nscf4c
                                                   | (Crystal, Basis, Eigvals)
                                                   +-> write_UNK_files -> tmp/UNK*
                                                   +-> write_input_files
                                                         +-> compute_pao_projectability
                                                         |     +-> UPF / Ggrid / ylmr / bessel / simpson
                                                         +-> save_scdm_provenance
                                                         +-> scdm.in + scdm.eig

scdm.in + scdm.eig + tmp/UNK* -> scdm_m -> scdm.amn + wan_r + Cubes
                                           |
               +---------------------------+--------------------------+
               |                           |                          |
 compute_scdm_only_bands          same function, Method=          run_wannier90
 Method='wi-scdm'                 'wi-scdm-f'                          |
               |                  transform -> project        .win/.eig/.amn
          project H(k)            -> FFT -> inverse               -> W90 -pp
               |                           |                     -> .nnkp
     centered-grid FFT             physical eigenvalues         -> write_w90_mmn
               |                           |                     -> W90
         interpolate                       |                     -> _hr.dat
               |                           |                          |
               +---------------------------+------------- compare_w90_bands
                                           |
                               calculate_err -> calculate_band_error
                                           |
                                  err.txt -> collect_eta_summary
```

Material launchers dispatch a validated stage name with `feval`. Only the selected case directory is placed on the path. `compare_w90_bands` and `run_eta_only` can use KSSOLV to regenerate missing DFT path references; this recovery path is an external dependency even during postprocessing. The Si MLWF comparison explicitly disables implicit recovery.

## Subspace construction

1. `write_input_files` accepts KSSOLV eigenvalues in Hartree, retaining all supplied bands, and writes absolute energies in eV. With an omitted target size, the raw PAO count sets `num_wann`; insufficient supplied bands raise an error.
2. PAOs are built from selected UPF `PP_PSWFC/PP_CHI.*` radial data and atomic positions. Orthogonalization uses an SVD rank test at each k. Count-only preflight reads UPF channels without evaluating wavefunction projectabilities. The optional fifth argument of `compute_pao_projectability` reads a compatible Basis directly instead of UNK text.
3. Manual target size plus valid manual `mu/sigma` skips PAO construction and fitting; isolated mode with an explicit target also skips them. Otherwise the necessary count and/or projectabilities are evaluated. Automatic erfc **and Gaussian** modes use the existing constrained erfc-fit policy; fitted `mu` is shifted by `-3*sigma_fit`. The Gaussian SCDM amplitude itself remains Gaussian. Nonfinite manual values or nonpositive manual sigma invoke the inherited warning/fallback policy (0/2 eV); supply valid parameters.
4. `scdm_m` chooses the input fractional k point minimizing `sum(abs(k),2)`; ties use input order. It does not perform a periodic nearest-Gamma search or an additional Gamma calculation.
5. The amplitude filter is `1` (isolated), `0.5*erfc((E-mu)/sigma)` (erfc), or `exp(-((E-mu)/sigma)^2)` (Gaussian). It multiplies the anchor wavefunctions directly before QRCP of their conjugate transpose.
6. Selected spin-space rows determine spatial centers `x_c = x - round(x)`. At each k, sampled cell-periodic states are multiplied by `exp(2*pi*i*x_c*k)` and the amplitude filter. Their conjugate transpose is factored by SVD, producing `A(k)=U*V'` with orthonormal columns. Real-space reconstruction restores the Bloch phase and averages over the uniform mesh.

For SOC, scalar PAOs are tensored with two spin components and the raw count doubles; this is not a separate j,mj recoupled atomic-orbital basis. Projectabilities use the identity plane-wave metric. PAW/ultrasoft overlap augmentation is not implemented. Norm-conserving UPF protocols obtain their data externally.

## Interpolation modes

| Mode | Entry | Operation and returned real-space matrix |
| --- | --- | --- |
| WI-SCDM | `compute_scdm_only_bands(...,'Method','wi-scdm')` | Lowdin `U=A*(A'*A)^(-1/2)`, `H_W=U'*diag(E)*U`, Fourier interpolation and diagonalization. `H_R` represents the physical projected Hamiltonian. |
| WI-SCDM-f | Same entry, `'Method','wi-scdm-f'` | Transform DFT energies, project/FFT/interpolate the transformed Hamiltonian, diagonalize, then invert the transform. `H_R` represents **f(H)**, not physical H. |
| SCDM+MLWF (`wlmf`) | `run_wannier90`, then `compare_w90_bands` | External Wannier90 initialized by SCDM. Default `dis_num_iter=0` preserves the supplied subspace while optimizing its gauge. HR interpolation uses degeneracies and optional WS-distance data. |

`read_scdm_uniform_kgrid` validates a complete Cartesian-product uniform mesh, its index mapping and mesh origin `delta`. Shifted meshes obey `k=m./N+delta` modulo reciprocal vectors. Returned coefficients for centered integer R contain `exp(-2*pi*i*delta.R)`. Evaluation uses origin-adjusted coordinates and cosine symmetrization of even-grid Nyquist boundaries. This is an integer FFT cell, not a metric-defined Wigner-Seitz construction. Two legacy Cartesian-coordinate conversions use cell metadata.

Transform defaults are `n=0.5`, `a=4*range(E_top)` (with a flat-band fallback), and `epsilon=max(E_top)` from the uniform-grid top band. Higher-level `compare_scdm_bands` instead supplies the largest finite DFT **path** energy for epsilon when unspecified. Retain this wrapper difference when comparing results. Record `transform_info` or set explicit parameters. No external Python `transform.py` is required.

The inverse uses at most 50 step-limited Newton iterations with transformed residual tolerance `1e-14`. The transform is flat above epsilon. Positive transformed eigenvalues are clipped to epsilon; diagnostics describe inversion, not material accuracy. Near the plateau, a small transformed residual need not imply a comparably small physical-energy error.

## Outputs and diagnostics

See [FILE_FORMATS.md](FILE_FORMATS.md) for shapes/references and [ERROR_METRICS_VITALE2020.md](ERROR_METRICS_VITALE2020.md) for metric sampling dependence.

Scalar Cubes contain the phase-aligned real part and can lose a significant imaginary component; the writer warns when this happens. Spinor Cubes contain total amplitude `sqrt(|up|^2+|down|^2)`, without spinor phase. MATLAB `wan_r` retains the complex two-component array. Cube spreads are real-space minimum-image estimates, distinct from Wannier90's reciprocal-space functional. Both direct methods share the same SCDM basis; there is no separately optimized SCDM-f orbital spread.

`save_scdm_provenance` records parameters, rank/fitting information and source/pseudopotential snapshots. It resolves functions at runtime with `which`, including external helpers. Generated archives contain external code and local paths and are excluded from the source package. `VALIDATION_RESULTS.json` records the current bundled-source static closure, supplemented by runtime file/class/command lookup review. Neither approach alone certifies every external calculation.
