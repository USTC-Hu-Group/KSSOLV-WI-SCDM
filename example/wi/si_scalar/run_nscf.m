function run_nscf
% Step 2. Uniform-grid eigenstates/UNK plus independent DFT path eigenvalues.
case_dir = si_case_setup;
diary(fullfile(case_dir, 'logs', 'run_nscf.log'));
cleanup = onCleanup(@() diary('off'));
s = load(fullfile(case_dir, 'scf_checkpoint.mat'), ...
    'cry', 'H', 'info', 'metadata', 'converged');
assert(s.converged, 'The SCF checkpoint is not converged.');
metadata = s.metadata;
rho = s.H.rho;
scf_info = s.info;
options_nscf = setksopt();
options_nscf.rho0 = rho;
options_nscf.X0 = [];
options_nscf.scfX = [];
options_nscf.iterscf = length(scf_info.Etotvec);
% autokpts resets fractional k points; finalize converts them exactly once.
crystal = set(s.cry, 'scfkpts', s.cry.kpts, ...
    'autokpts', metadata.nscf_autokpts, 'nbnd', metadata.nbnd_nscf);
crystal = finalize(crystal);
crystal.ppvar = []; % Rebuild nonlocal lookup tables for this k-point grid.
rng(metadata.random_seed, 'twister');
[crystal, BH, BX, info] = nscf4c(crystal, options_nscf);
assert(numel(info.Eigvals) == metadata.nbnd_nscf*512, ...
    'Uniform NSCF must return 512 k points and 12 bands.');
save(fullfile(case_dir, 'nscf_checkpoint.mat'), 'crystal', 'BH', ...
    'BX', 'info', 'options_nscf', 'rho', 'metadata', '-v7.3');
write_UNK_files(BX, crystal);

[~, k_path_fractional, path_distance_bohr_inv] = si_band_path(crystal.supercell);
crystal_path = set(crystal, 'kpts', k_path_fractional, ...
    'scfkpts', crystal.kpts);
crystal_path = finalize(crystal_path);
crystal_path.ppvar = [];
[crystal_path, ~, ~, info_path] = nscf4c(crystal_path, options_nscf);
energyBands = reshape(info_path.Eigvals, metadata.nbnd_nscf, []);
assert(size(energyBands, 2) == size(k_path_fractional, 1), ...
    'DFT path eigenvalues do not match the interpolation path.');
save(fullfile(case_dir, 'nscf_path_checkpoint.mat'), 'crystal_path', ...
    'energyBands', 'metadata', 'k_path_fractional', ...
    'path_distance_bohr_inv', '-v7.3');
energy_relative_ev = (energyBands.' - crystal.efermi)*27.211386246;
names = arrayfun(@(b) sprintf('DFT_band_%02d_eV', b), ...
    1:size(energy_relative_ev, 2), 'UniformOutput', false);
data = array2table([k_path_fractional, path_distance_bohr_inv, energy_relative_ev], ...
    'VariableNames', [{'k1_frac','k2_frac','k3_frac','distance_bohr_inv'}, names]);
writetable(data, fullfile(case_dir, 'results', 'dft_bands.csv'));
fprintf('=== Si NSCF complete: 512 mesh points, %d path points ===\n', ...
    size(k_path_fractional, 1));
end
