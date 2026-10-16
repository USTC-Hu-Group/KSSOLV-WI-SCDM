function result = run_wlmf(varargin)
%RUN_WLMF Run Wannier90 localization and compare MLWF/Wannier bands to DFT.
%
% The method label written to err.txt is intentionally "wlmf" to match the
% current benchmark naming in this folder. It corresponds to the Wannier90
% localized/MLWF Hamiltonian built from SCDM projections.

ensure_dir('logs');
ensure_dir('results');
ensure_dir('wannier90_post');
diary(fullfile('logs', 'run_wlmf.log'));
cleanup = onCleanup(@() diary('off'));
validate_vitale_error_runtime;
enable_soc_blochham_compat;

metadata = struct();
if isfile('nscf_checkpoint.mat')
    loaded = load('nscf_checkpoint.mat', 'metadata');
    if isfield(loaded, 'metadata')
        metadata = loaded.metadata;
    end
end
case_name = get_metadata_text(metadata, 'case_name', get_current_folder_name());
fprintf('=== WLMF/Wannier90 start: %s ===\n', case_name);

w90_command = get_env_text('W90_COMMAND', 'wannier90.x');
num_iter = get_env_number('W90_NUM_ITER', 10000);
dis_num_iter = get_env_number('W90_DIS_NUM_ITER', 0);
conv_tol = get_env_optional_number('W90_CONV_TOL');
conv_window = get_env_number('W90_CONV_WINDOW', 3);
kmesh_tol = get_env_number('W90_KMESH_TOL', 1e-8);
do_compare = get_env_logical('W90_COMPARE_BANDS', true);
plot_band_error = get_env_logical('W90_PLOT_BAND_ERROR', false);
band_error_mode = get_metadata_text(metadata, 'band_error_mode', 'isolated');

result = run_wannier90( ...
    'WorkingDirectory', pwd, ...
    'SeedName', 'wannier90', ...
    'OutputDirectory', 'wannier90_post', ...
    'W90Command', w90_command, ...
    'NumIter', num_iter, ...
    'DisNumIter', dis_num_iter, ...
    'ConvTol', conv_tol, ...
    'ConvWindow', conv_window, ...
    'KmeshTol', kmesh_tol, ...
    'ErrorOnFailure', true, ...
    'WannierPlot', true, ...
    'WannierPlotSupercell', get_env_number('W90_PLOT_SUPERCELL', 1), ...
    'WannierPlotSpinorMode', 'total', ...
    'Force', true);

if do_compare
    ensure_path_band_checkpoint(metadata);
    try
        result.band_compare = compare_w90_bands( ...
            'WorkingDirectory', pwd, ...
            'SeedName', 'wannier90', ...
            'W90Directory', 'wannier90_post', ...
            'Method', 'wlmf', ...
            'BandErrorMode', band_error_mode, ...
            'PlotBandError', plot_band_error, ...
            'ErrorPlotPrefix', fullfile('results', 'wlmf_band_error_metric'));
        if ~has_vitale_band_error('err.txt', 'wlmf', band_error_mode, 1, 0.1)
            error('run_wlmf:MissingVitaleMetrics', ...
                'Missing current wlmf eta/eta_max record. Check the active branch src/WI path.');
        end
    catch compare_error
        warning('run_wlmf:CompareFailed', ...
            'Wannier90/WLMF band comparison failed: %s', compare_error.message);
        result.band_compare_error = compare_error;
    end
end

save('wlmf_checkpoint.mat', 'result', '-v7.3');
if isfield(result, 'band_compare_error')
    rethrow(result.band_compare_error);
end
fprintf('=== WLMF/Wannier90 done: %s ===\n', case_name);
end

function ensure_path_band_checkpoint(metadata)
if isfile('nscf_path_checkpoint.mat') || ~isfile('nscf_checkpoint.mat')
    return;
end

try
    nscf_data = load('nscf_checkpoint.mat', ...
        'crystal', 'metadata', 'options_nscf', 'rho');
    if isfield(nscf_data, 'metadata')
        metadata = nscf_data.metadata;
    end
    if ~has_path_kpoints(metadata) || ~isfield(nscf_data, 'options_nscf') || ...
            ~isfield(nscf_data, 'rho')
        return;
    end

    [crystal_path, energyBands] = run_path_bands_with_refreshed_ppvar( ...
        nscf_data.crystal, nscf_data.options_nscf, ...
        metadata.path_kpoints, nscf_data.rho);
    save('nscf_path_checkpoint.mat', 'crystal_path', ...
        'energyBands', 'metadata', '-v7.3');
    fprintf('Recovered missing path-band checkpoint: %s\n', ...
        fullfile(pwd, 'nscf_path_checkpoint.mat'));
catch path_error
    warning('run_wlmf:PathBandsFailed', ...
        'Path-band NSCF recovery failed: %s', path_error.message);
end
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

function tf = has_path_kpoints(metadata)
tf = isstruct(metadata) && isfield(metadata, 'path_kpoints') && ...
    iscell(metadata.path_kpoints) && size(metadata.path_kpoints, 1) >= 2 && ...
    size(metadata.path_kpoints, 2) >= 5;
end

function value = get_env_text(name, default_value)
value = getenv(name);
if isempty(strtrim(value))
    value = default_value;
end
end

function value = get_env_number(name, default_value)
text = getenv(name);
if isempty(strtrim(text))
    value = default_value;
    return;
end
value = str2double(text);
if ~isfinite(value)
    warning('run_wlmf:BadEnvironmentNumber', ...
        'Environment variable %s=%s is not numeric; using %g.', ...
        name, text, default_value);
    value = default_value;
end
end

function value = get_env_optional_number(name)
text = getenv(name);
if isempty(strtrim(text)) || strcmpi(strtrim(text), 'auto')
    value = [];
    return;
end
value = str2double(text);
if ~isfinite(value)
    warning('run_wlmf:BadEnvironmentNumber', ...
        'Environment variable %s=%s is not numeric; using auto.', name, text);
    value = [];
end
end

function tf = get_env_logical(name, default_value)
text = getenv(name);
if isempty(strtrim(text))
    tf = default_value;
else
    tf = any(strcmpi(strtrim(text), {'1', 'true', 'yes', 'on'}));
end
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

function name = get_current_folder_name()
[~, name] = fileparts(pwd);
end

function ensure_dir(dirname)
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
end
