# Band-error definitions and records

The retained implementation labels its formulas `vitale2020-v1`, referring to Vitale et al., *Automated high-throughput Wannierisation*, npj Computational Materials 6, 66 (2020), [DOI: 10.1038/s41524-020-0312-y](https://doi.org/10.1038/s41524-020-0312-y). This guide describes the actual software convention; it does not certify equivalence to an independently obtained implementation of that article.

`calculate_band_error` takes `[band,kpoint]` arrays in eV relative to the same DFT Fermi energy. Bands/k points must correspond; it does not perform band assignment or fit an energy offset.

For isolated bands, the implementation follows the unnormalized printed expression:

```text
delta = E_DFT - E_Wan
eta = sqrt(sum(delta.^2))
eta_max = max(abs(delta))
```

The sum covers all supplied bands and k points. Thus isolated eta is sampling-dependent and is **not RMS**. Ordinary RMS is reported separately. Do not reuse an RMS threshold for this eta or pool values with different sample sizes as if normalized.

For entangled bands:

```text
f(E) = 1 / (1 + exp((E-Nu)/Tau))
w = sqrt(f(E_DFT) * f(E_Wan))
eta = sqrt(sum(w .* delta.^2) / sum(w))
eta_max = max(w .* abs(delta))
```

`eta_max` uses w, not sqrt(w), and differs from the unweighted maximum error. Weights depend on both energy sets. Log-space Fermi-Dirac evaluation avoids premature underflow; an all-zero-weight comparison raises an error. Defaults are `Nu=1.0` eV relative to the common DFT Fermi energy and `Tau=0.1` eV. `Nu/Tau` are independent of the SCDM construction filter's absolute-energy `mu/sigma`.

```matlab
[eta, detail] = calculate_band_error(E_DFT, E_Wan, 'entangled', ...
    'Nu', 1, 'Tau', 0.1);
% bands: [kpoint,band], relative eV; energyBands: [band,kpoint], Hartree.
[mae, maxerr, rms, detail] = calculate_err(bands, energyBands, crystal, ...
    'BandErrorMode', 'entangled', 'Method', 'wi-scdm');
```

`calculate_err` converts the KSSOLV reference with `crystal.efermi`, writes versioned records to `err.txt`, and retains MAE, RMS and unweighted MAX in addition to eta/eta_max. `ReferenceBandStart` and `DiscardTopBands` explicitly control selection; default evaluation does not discard top bands. Supplying obsolete metric `Mu/Sigma` options is rejected rather than silently reusing construction parameters.

`has_vitale_band_error` and `run_eta_only` require a complete matching record, including formula version, exact method, mode and `Nu/Tau`, before skipping work. Existing legacy metrics are not relabeled as current results. `run_eta_only('WorkingDirectory',case_dir)` can recompute metrics from saved data; if a DFT path reference is missing, its recovery path can invoke KSSOLV NSCF.

`collect_eta_summary(case_root)` summarizes current records with the default evaluation window and pairs comparable method results. It retains eta/eta_max and MAE/RMS/MAX, and rejects pairing with different mode, reference selection or comparison dimensions. Pairwise delta eta is MLWF minus the respective direct method; positive values mean a larger MLWF error under this metric. No external Python batch collector or Slurm script is needed by this release.

Tests cover analytic weights, the isolated eta/RMS distinction, weighted maximum, exchange symmetry, a common energy shift, numerical tails, Hartree conversion, filter/metric independence, explicit reference trimming and cache/summary selection. See [tests/README.md](../../tests/wi/README.md) for commands and validation boundaries.
