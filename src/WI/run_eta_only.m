function result = run_eta_only(varargin)
%RUN_ETA_ONLY Recompute eta reports for an existing SCDM/WLMF case.
%
% Reuses SCF/NSCF/SCDM/Wannier90 outputs and appends Vitale 2020 metrics.
% If reference path bands are missing, the existing path-NSCF recovery is
% used. Nu (default 1 eV above DFT E_F) and Tau (0.1 eV) are independent of
% scdm.in. Legacy eta records never satisfy the completion check.

parser = inputParser;
parser.FunctionName = mfilename;
parser.PartialMatching = false;
addParameter(parser, 'WorkingDirectory', pwd, @(x) ischar(x) || isstring(x));
addParameter(parser, 'ComputeSCDM', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'ComputeSCDMF', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'ComputeWLMF', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'Force', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'PlotBandError', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'Nu', 1.0, @is_scalar_real);
addParameter(parser, 'Tau', 0.1, @(x) is_scalar_real(x) && x > 0);
parse(parser, varargin{:});

case_dir = char(string(parser.Results.WorkingDirectory));
compute_scdm = logical(parser.Results.ComputeSCDM);
compute_scdm_f = logical(parser.Results.ComputeSCDMF);
compute_wlmf = logical(parser.Results.ComputeWLMF);
force = logical(parser.Results.Force);
plot_band_error = logical(parser.Results.PlotBandError);
nu = double(parser.Results.Nu);
tau = double(parser.Results.Tau);

if ~isfolder(case_dir)
    error('run_eta_only:MissingCaseDirectory', ...
        'WorkingDirectory does not exist: %s', case_dir);
end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(case_dir);

ensure_dir('logs');
ensure_dir('results');
diary_file = fullfile('logs', 'run_eta_only.log');
diary(diary_file);
cleanup_diary = onCleanup(@() diary('off'));

fprintf('=== ETA-only recomputation start: %s ===\n', case_dir);
fprintf('Started at %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

nscf_file = fullfile(case_dir, 'nscf_checkpoint.mat');
if ~isfile(nscf_file)
    error('run_eta_only:MissingNSCFCheckpoint', ...
        'Missing nscf_checkpoint.mat in %s.', case_dir);
end

nscf_data = load(nscf_file, 'crystal', 'metadata', 'options_nscf', 'rho');
if ~isfield(nscf_data, 'crystal')
    error('run_eta_only:MissingCrystal', ...
        'nscf_checkpoint.mat does not contain variable crystal.');
end
crystal = nscf_data.crystal;
metadata = get_loaded_metadata(nscf_data);
band_error_mode = validatestring( ...
    get_metadata_text(metadata, 'band_error_mode', 'isolated'), ...
    {'isolated', 'entangled'}, mfilename, 'band_error_mode');

result = struct();
result.case_dir = case_dir;
result.band_error_mode = band_error_mode;
result.metric_definition = 'vitale2020-v1';
result.nu_relative_fermi_eV = nu;
result.tau_eV = tau;
result.scdm = [];
result.scdm_f = [];
result.wlmf = [];
result.scdm_error = [];
result.scdm_f_error = [];
result.wlmf_error = [];

try
    ensure_path_band_checkpoint(nscf_data, metadata);
catch err
    warning('run_eta_only:PathBandCheckpointFailed', ...
        'Could not ensure nscf_path_checkpoint.mat: %s', err.message);
end

if ~isfile('nscf_path_checkpoint.mat')
    error('run_eta_only:MissingPathBands', ...
        'Missing nscf_path_checkpoint.mat; eta cannot be recomputed.');
end
path_data = load('nscf_path_checkpoint.mat', 'energyBands');
if ~isfield(path_data, 'energyBands')
    error('run_eta_only:MissingEnergyBands', ...
        'nscf_path_checkpoint.mat does not contain energyBands.');
end
energyBands = path_data.energyBands;

if compute_scdm
    if ~force && has_vitale_band_error('err.txt', 'wi-scdm', band_error_mode, nu, tau)
        fprintf('Skipping wi-scdm eta: err.txt already contains it.\n');
    elseif has_scdm_prerequisites()
        try
            fprintf('Recomputing wi-scdm eta from existing scdm.amn/scdm.eig/scdm.in...\n');
            [bands_scdm, H_R] = compare_scdm_bands( ...
                metadata.path_kpoints, energyBands, crystal, 'scdm.amn', ...
                'Method', 'wi-scdm');
            [err_aver, err_max, err_rms, detail] = calculate_err( ...
                bands_scdm, energyBands, crystal, ...
                'Method', 'wi-scdm', ...
                'BandErrorMode', band_error_mode, ...
                'Nu', nu, 'Tau', tau, ...
                'PlotBandError', plot_band_error, ...
                'ErrorPlotPrefix', fullfile('results', 'eta_only_wi_scdm_band_error_metric'));
            result.scdm = struct( ...
                'bands', bands_scdm, ...
                'H_R', H_R, ...
                'err_aver', err_aver, ...
                'err_max', err_max, ...
                'err_rms', err_rms, ...
                'detail', detail);
        catch err
            warning('run_eta_only:SCDMError', ...
                'wi-scdm eta recomputation failed: %s', err.message);
            result.scdm_error = err;
        end
    else
        warning('run_eta_only:MissingSCDMFiles', ...
            'Skipping wi-scdm eta because scdm.amn/scdm.eig/scdm.in is incomplete.');
    end
end

if compute_scdm_f
    if ~force && has_vitale_band_error('err.txt', 'wi-scdm-f', band_error_mode, nu, tau)
        fprintf('Skipping wi-scdm-f eta: err.txt already contains it.\n');
    elseif has_scdm_prerequisites()
        try
            fprintf('Recomputing wi-scdm-f eta from existing scdm.amn/scdm.eig/scdm.in...\n');
            [bands_scdm_f, H_R_f] = compare_scdm_bands( ...
                metadata.path_kpoints, energyBands, crystal, 'scdm.amn', ...
                'Method', 'wi-scdm-f');
            [err_aver_f, err_max_f, err_rms_f, detail_f] = calculate_err( ...
                bands_scdm_f, energyBands, crystal, ...
                'Method', 'wi-scdm-f', ...
                'BandErrorMode', band_error_mode, ...
                'Nu', nu, 'Tau', tau, ...
                'PlotBandError', plot_band_error, ...
                'ErrorPlotPrefix', fullfile('results', 'eta_only_wi_scdm_f_band_error_metric'));
            result.scdm_f = struct( ...
                'bands', bands_scdm_f, ...
                'H_R', H_R_f, ...
                'err_aver', err_aver_f, ...
                'err_max', err_max_f, ...
                'err_rms', err_rms_f, ...
                'detail', detail_f);
        catch err
            warning('run_eta_only:SCDMFError', ...
                'wi-scdm-f eta recomputation failed: %s', err.message);
            result.scdm_f_error = err;
        end
    else
        warning('run_eta_only:MissingSCDMFiles', ...
            'Skipping wi-scdm-f eta because scdm.amn/scdm.eig/scdm.in is incomplete.');
    end
end

if compute_wlmf
    if ~force && has_vitale_band_error('err.txt', 'wlmf', band_error_mode, nu, tau)
        fprintf('Skipping wlmf eta: err.txt already contains it.\n');
    elseif has_wlmf_prerequisites()
        try
            fprintf('Recomputing wlmf eta from existing wannier90_post/wannier90_hr.dat...\n');
            result.wlmf = compare_w90_bands( ...
                'WorkingDirectory', case_dir, ...
                'SeedName', 'wannier90', ...
                'W90Directory', 'wannier90_post', ...
                'Method', 'wlmf', ...
                'BandErrorMode', band_error_mode, ...
                'Nu', nu, 'Tau', tau, ...
                'PlotBandError', plot_band_error, ...
                'ErrorPlotPrefix', fullfile('results', 'eta_only_wlmf_band_error_metric'));
        catch err
            warning('run_eta_only:WLMFError', ...
                'wlmf eta recomputation failed: %s', err.message);
            result.wlmf_error = err;
        end
    else
        warning('run_eta_only:MissingWLMFFiles', ...
            'Skipping wlmf eta because wannier90_post/wannier90_hr.dat is missing.');
    end
end

save('eta_only_checkpoint.mat', 'result', '-v7.3');
fprintf('ETA-only checkpoint saved: %s\n', fullfile(case_dir, 'eta_only_checkpoint.mat'));
fprintf('Finished at %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf('=== ETA-only recomputation done: %s ===\n', case_dir);

clear cleanup_diary cleanup_dir;
end

function metadata = get_loaded_metadata(nscf_data)
if isfield(nscf_data, 'metadata') && isstruct(nscf_data.metadata)
    metadata = nscf_data.metadata;
else
    metadata = struct();
end
end

function ensure_path_band_checkpoint(nscf_data, metadata)
if isfile('nscf_path_checkpoint.mat')
    return;
end
if ~has_path_kpoints(metadata) || ~isfield(nscf_data, 'options_nscf') || ...
        ~isfield(nscf_data, 'rho')
    return;
end

fprintf('nscf_path_checkpoint.mat is missing; recomputing reference path bands.\n');
[crystal_path, energyBands] = run_path_bands_with_refreshed_ppvar( ...
    nscf_data.crystal, nscf_data.options_nscf, metadata.path_kpoints, ...
    nscf_data.rho);
save('nscf_path_checkpoint.mat', 'crystal_path', 'energyBands', ...
    'metadata', '-v7.3');
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

function tf = is_scalar_real(value)
tf = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end

function tf = has_scdm_prerequisites()
tf = isfile('scdm.amn') && isfile('scdm.eig') && isfile('scdm.in');
end

function tf = has_wlmf_prerequisites()
tf = isfile(fullfile('wannier90_post', 'wannier90_hr.dat'));
end

function tf = has_path_kpoints(metadata)
tf = isstruct(metadata) && isfield(metadata, 'path_kpoints') && ...
    iscell(metadata.path_kpoints) && size(metadata.path_kpoints, 2) >= 5 && ...
    size(metadata.path_kpoints, 1) >= 2;
end

function value = get_metadata_text(metadata, field_name, default_value)
if isstruct(metadata) && isfield(metadata, field_name) && ...
        ~isempty(metadata.(field_name))
    value = char(string(metadata.(field_name)));
else
    value = default_value;
end
end

function ensure_dir(dirname)
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
end
