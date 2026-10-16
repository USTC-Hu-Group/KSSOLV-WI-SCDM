function run_scdm
ensure_dir('logs');
ensure_dir('results');
diary(fullfile('logs', 'run_scdm.log'));
cleanup = onCleanup(@() diary('off'));
validate_vitale_error_runtime;

enable_soc_blochham_compat;
kssolvpptype('default');
load('nscf_checkpoint.mat', 'crystal', 'info', 'metadata', ...
    'options_nscf', 'rho');
fprintf('=== SCDM start: %s ===\n', metadata.case_name);

band_error_mode = get_metadata_text(metadata, 'band_error_mode', 'isolated');
write_input_files(crystal, info.Eigvals, [], metadata.focc);

[wan_r, Amn, cube_files] = scdm_m('scdm.in', 'scdm.amn');
save('scdm_checkpoint.mat', 'wan_r', 'Amn', 'cube_files', ...
    'metadata', '-v7.3');

if should_compute_spread()
    try
        [spread_table, spread_summary] = calculate_scdm_spread( ...
            'CubePattern', 'scdm_*.cube', ...
            'OutputPrefix', fullfile('results', [metadata.case_name '_scdm_spread']));
        save('spread_checkpoint.mat', 'spread_table', 'spread_summary', ...
            'metadata', '-v7.3');
    catch spread_error
        warning('run_scdm:SpreadFailed', ...
            'SCDM spread calculation failed: %s', spread_error.message);
    end
end

if ~isfile('nscf_path_checkpoint.mat') && ...
        isfield(metadata, 'run_path_bands') && metadata.run_path_bands
    try
        [crystal_path, energyBands] = run_path_bands_with_refreshed_ppvar( ...
            crystal, options_nscf, metadata.path_kpoints, rho);
        save('nscf_path_checkpoint.mat', 'crystal_path', ...
            'energyBands', 'metadata', '-v7.3');
        fprintf('Recovered missing path-band checkpoint: %s\n', ...
            fullfile(pwd, 'nscf_path_checkpoint.mat'));
    catch path_error
        warning('run_scdm:PathBandsFailed', ...
            'Path-band NSCF recovery failed: %s', path_error.message);
    end
end

if isfile('nscf_path_checkpoint.mat')
    load('nscf_path_checkpoint.mat', 'energyBands');

    [bands_scdm, H_R] = compare_scdm_bands( ...
        metadata.path_kpoints, energyBands, crystal, 'scdm.amn', ...
        'Method', 'wi-scdm');
    [err_aver, err_max, err_rms, band_error_detail] = calculate_err( ...
        bands_scdm, energyBands, crystal, ...
        'Method', 'wi-scdm', ...
        'BandErrorMode', band_error_mode, ...
        'ErrorPlotPrefix', fullfile('results', 'wi_scdm_band_error_metric'));

    [bands_scdm_f, H_R_f] = compare_scdm_bands( ...
        metadata.path_kpoints, energyBands, crystal, 'scdm.amn', ...
        'Method', 'wi-scdm-f');
    [err_aver_f, err_max_f, err_rms_f, band_error_detail_f] = calculate_err( ...
        bands_scdm_f, energyBands, crystal, ...
        'Method', 'wi-scdm-f', ...
        'BandErrorMode', band_error_mode, ...
        'ErrorPlotPrefix', fullfile('results', 'wi_scdm_f_band_error_metric'));

    save('band_compare_checkpoint.mat', 'bands_scdm', 'H_R', ...
        'bands_scdm_f', 'H_R_f', 'energyBands', ...
        'err_aver', 'err_max', 'err_rms', 'band_error_detail', ...
        'err_aver_f', 'err_max_f', 'err_rms_f', 'band_error_detail_f', ...
        '-v7.3');
    for method = {'wi-scdm', 'wi-scdm-f'}
        if ~has_vitale_band_error('err.txt', method{1}, band_error_mode, 1, 0.1)
            error('run_scdm:MissingVitaleMetrics', ...
                'Missing current eta/eta_max record for %s. Check the active branch src/WI path.', method{1});
        end
    end
else
    error('run_scdm:MissingPathBands', ...
        'nscf_path_checkpoint.mat is missing; eta/eta_max require DFT path bands.');
end

fprintf('=== SCDM done: %s ===\n', metadata.case_name);
end

function tf = should_compute_spread()
value = getenv('SCDM_COMPUTE_SPREAD');
tf = strcmp(value, '1') || strcmpi(value, 'true');
end

function [crystal_path, energyBands] = run_path_bands_with_refreshed_ppvar( ...
        crystal, options_nscf, path_kpoints, rho)
symkpts = reshape([path_kpoints{:, 1:3}], [], 3);
nkpts = [path_kpoints{:, 4}];
kpoints = [path_kpoints{1, 1:3}];
for i = 1:size(path_kpoints, 1) - 1
    nk = nkpts(i) - 1;
    dk = (symkpts(i + 1, 1:3) - symkpts(i, 1:3)) .* ((1:nk)' / nk);
    kpoints = [kpoints; symkpts(i, 1:3) + dk]; %#ok<AGROW>
end

options_nscf.rho0 = rho;
options_nscf.scfX = [];
options_nscf.X0 = [];
crystal_path = set(crystal, 'kpts', kpoints, 'scfkpts', crystal.kpts);
crystal_path = finalize(crystal_path);
crystal_path.ppvar = [];

[crystal_path, ~, ~, info_path] = nscf4c(crystal_path, options_nscf);
energyBands = reshape(info_path.Eigvals, crystal_path.nbnd, []);
end

function value = get_metadata_text(metadata, field_name, default_value)
value = default_value;
if isstruct(metadata) && isfield(metadata, field_name)
    candidate = metadata.(field_name);
    if ischar(candidate) || (isstring(candidate) && isscalar(candidate))
        value = char(string(candidate));
    end
end
end

function ensure_dir(dirname)
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
end
