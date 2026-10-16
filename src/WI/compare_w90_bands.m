function result = compare_w90_bands(varargin)
%COMPARE_W90_BANDS Compare Wannier90 HR interpolation with KSSOLV DFT bands.
%
% result = compare_w90_bands('WorkingDirectory', case_dir)
%
% The function reads wannier90_post/wannier90_hr.dat, interpolates it on the
% exact k-path stored in metadata.path_kpoints, and compares against
% nscf_path_checkpoint.mat. Energies are reported in eV relative to E_F.

parser = inputParser;
parser.FunctionName = mfilename;
parser.PartialMatching = false;
addParameter(parser, 'WorkingDirectory', pwd, @(x) ischar(x) || isstring(x));
addParameter(parser, 'W90Directory', 'wannier90_post', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SeedName', 'wannier90', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Method', 'wannier90', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(parser, 'OutputPrefix', fullfile('results', 'compare_wannier90_bands'), ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'Visible', 'off', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ReferenceBandStart', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && ...
    x >= 1 && x == round(x));
addParameter(parser, 'RecomputeReferenceIfMissing', true, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'BandErrorMode', 'isolated', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(parser, 'Mu', [], @is_empty_or_scalar_real);
addParameter(parser, 'Sigma', [], @is_empty_or_scalar_real);
addParameter(parser, 'Nu', 1.0, @(x) ~isempty(x) && is_empty_or_scalar_real(x));
addParameter(parser, 'Tau', 0.1, ...
    @(x) ~isempty(x) && is_empty_or_scalar_real(x) && x > 0);
addParameter(parser, 'PlotBandError', false, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'plot_band_error', [], ...
    @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(parser, 'ErrorPlotPrefix', fullfile('results', 'wannier90_band_error_metric'), ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(parser, varargin{:});
if ~isempty(parser.Results.Mu) || ~isempty(parser.Results.Sigma)
    error('compare_w90_bands:SCDMParametersNotEvaluationParameters', ...
        'Remove Mu/Sigma from error evaluation; use independent Nu and Tau.');
end

case_dir = char(parser.Results.WorkingDirectory);
w90_dir = char(parser.Results.W90Directory);
seedname = char(parser.Results.SeedName);
method_label = char(string(parser.Results.Method));
output_prefix = char(parser.Results.OutputPrefix);
visible = char(parser.Results.Visible);
reference_band_start = double(parser.Results.ReferenceBandStart);
recompute_reference = logical(parser.Results.RecomputeReferenceIfMissing);
band_error_mode = validatestring( ...
    lower(char(string(parser.Results.BandErrorMode))), ...
    {'isolated', 'entangled'}, mfilename, 'BandErrorMode');
plot_band_error = logical(parser.Results.PlotBandError);
if ~isempty(parser.Results.plot_band_error)
    plot_band_error = logical(parser.Results.plot_band_error);
end
error_plot_prefix = char(string(parser.Results.ErrorPlotPrefix));

if ~isfolder(case_dir)
    error('compare_w90_bands:MissingWorkingDirectory', ...
        'WorkingDirectory does not exist: %s', case_dir);
end
if ~is_absolute_path(w90_dir)
    w90_dir = fullfile(case_dir, w90_dir);
end
if ~is_absolute_path(output_prefix)
    output_prefix = fullfile(case_dir, output_prefix);
end
output_dir = fileparts(output_prefix);
if ~isempty(output_dir) && ~isfolder(output_dir)
    mkdir(output_dir);
end

nscf_file = fullfile(case_dir, 'nscf_checkpoint.mat');
if ~isfile(nscf_file)
    error('compare_w90_bands:MissingNSCFCheckpoint', ...
        'Missing nscf_checkpoint.mat in %s.', case_dir);
end
nscf_data = load(nscf_file, 'crystal', 'metadata', 'options_nscf', 'rho');
crystal = nscf_data.crystal;
metadata = nscf_data.metadata;

if ~has_path_kpoints(metadata)
    error('compare_w90_bands:MissingKPath', ...
        'metadata.path_kpoints is missing or malformed in %s.', nscf_file);
end

reference_file = fullfile(case_dir, 'nscf_path_checkpoint.mat');
if isfile(reference_file)
    reference_data = load(reference_file, 'energyBands');
    energy_bands_hartree = reference_data.energyBands;
elseif recompute_reference && isfield(nscf_data, 'options_nscf') && ...
        isfield(nscf_data, 'rho')
    warning('compare_w90_bands:RecomputeReference', ...
        'nscf_path_checkpoint.mat is missing; recomputing path bands.');
    [crystal_path, energyBands] = recompute_reference_bands( ...
        crystal, nscf_data.options_nscf, metadata.path_kpoints, ...
        nscf_data.rho);
    energy_bands_hartree = energyBands;
    save(reference_file, 'crystal_path', 'energyBands', 'metadata', '-v7.3');
else
    error('compare_w90_bands:MissingReferenceBands', ...
        ['Missing nscf_path_checkpoint.mat, and the reference bands cannot ', ...
         'be recomputed from nscf_checkpoint.mat.']);
end

hr_file = fullfile(w90_dir, [seedname, '_hr.dat']);
if ~isfile(hr_file)
    error('compare_w90_bands:MissingHR', ...
        'Missing Wannier90 HR file: %s. Run run_wannier90 first.', hr_file);
end

[k_path, path_dist, node_locs, path_labels] = ...
    build_k_path(metadata.path_kpoints);
hr = read_w90_hr(hr_file);
bands_w90_abs = interpolate_w90_hr(hr, k_path);
efermi_ev = get_fermi_energy(crystal) * 27.211386246;
bands_w90 = bands_w90_abs - efermi_ev;

[energy_ref, reference_band_indices] = orient_and_select_reference( ...
    energy_bands_hartree, size(bands_w90), reference_band_start);
energy_ref = (energy_ref - get_fermi_energy(crystal)) * 27.211386246;

old_dir = pwd;
cleanup_obj = onCleanup(@() cd(old_dir));
cd(case_dir);
[err_aver, err_max, err_rms, band_error_detail] = calculate_err( ...
    bands_w90, energy_bands_hartree, crystal, ...
    'Method', method_label, ...
    'ReferenceBandStart', reference_band_start, ...
    'BandErrorMode', band_error_mode, ...
    'Nu', parser.Results.Nu, ...
    'Tau', parser.Results.Tau, ...
    'PlotBandError', plot_band_error, ...
    'ErrorPlotPrefix', error_plot_prefix, ...
    'Visible', visible);

fig = figure('Color', 'w', 'Visible', visible, ...
    'Name', 'Wannier90 Band Comparison');
hold on;
plot(path_dist, bands_w90, 'b--', 'LineWidth', 1.4);
plot(path_dist, energy_ref, 'k-', 'LineWidth', 1.2);
xlim([0, max(path_dist)]);
all_values = [bands_w90(:); energy_ref(:)];
all_values = all_values(isfinite(all_values));
if ~isempty(all_values)
    ymin = min(all_values);
    ymax = max(all_values);
    ypad = max(0.05 * (ymax - ymin), 0.1);
    ylim([ymin - ypad, ymax + ypad]);
end
set(gca, 'XTick', node_locs);
set(gca, 'XTickLabel', path_labels);
for loc = node_locs
    xline(loc, 'k:', 'HandleVisibility', 'off');
end
grid on;
box on;
ylabel('Energy - E_F (eV)');
title('Wannier90 Interpolated Bands Compared with KSSOLV DFT');
legend({'Wannier90', 'KSSOLV DFT'}, 'Location', 'best');
hold off;

png_file = [output_prefix, '.png'];
fig_file = [output_prefix, '.fig'];
csv_file = [output_prefix, '.csv'];
mat_file = [output_prefix, '.mat'];
exportgraphics(fig, png_file, 'Resolution', 300);
savefig(fig, fig_file);
write_band_csv(csv_file, path_dist, bands_w90, energy_ref);

result = struct();
result.case_dir = case_dir;
result.w90_dir = w90_dir;
result.seedname = seedname;
result.method = method_label;
result.hr_file = hr_file;
result.png_file = png_file;
result.fig_file = fig_file;
result.csv_file = csv_file;
result.err_aver = err_aver;
result.err_max = err_max;
result.err_rms = err_rms;
result.band_error_mode = band_error_mode;
result.band_eta = band_error_detail.eta;
result.band_eta_max = band_error_detail.eta_max;
result.band_error_detail = band_error_detail;
result.reference_band_indices = reference_band_indices;
result.bands_w90 = bands_w90;
result.energy_ref = energy_ref;
result.path_dist = path_dist;
result.k_path = k_path;
save(mat_file, 'result', '-v7.3');

fprintf('Wannier90 band comparison saved to %s\n', png_file);
fprintf('Wannier90 band data saved to %s\n', csv_file);
end

function tf = is_empty_or_scalar_real(value)
tf = isempty(value) || ...
    (isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value));
end

function [crystal_path, energyBands] = recompute_reference_bands( ...
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

% Changing from the uniform NSCF mesh to a high-symmetry path invalidates the
% nonlocal pseudopotential interpolation table stored in ppvar. Rebuilding it
% avoids getvnl indexing past the table near path-boundary k points.
crystal_path.ppvar = [];

[crystal_path, ~, ~, info_path] = nscf4c(crystal_path, options_nscf);
energyBands = reshape(info_path.Eigvals, crystal_path.nbnd, []);
end

function tf = is_absolute_path(path_value)
path_value = char(path_value);
tf = startsWith(path_value, filesep) || ...
    ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once'));
end

function tf = has_path_kpoints(metadata)
tf = isstruct(metadata) && isfield(metadata, 'path_kpoints') && ...
    iscell(metadata.path_kpoints) && size(metadata.path_kpoints, 1) >= 2 && ...
    size(metadata.path_kpoints, 2) >= 5;
end

function [k_path, path_dist, node_locs, path_labels] = build_k_path(k_points)
num_nodes = size(k_points, 1);
path_nodes = cell2mat(k_points(:, 1:3));
path_labels = k_points(:, 5);
segment_counts = cell2mat(k_points(1:end-1, 4));
if any(segment_counts < 2) || any(segment_counts ~= round(segment_counts))
    error('compare_w90_bands:BadKPath', ...
        'Each k-path segment must contain an integer number of at least 2 points.');
end

nk_path = sum(segment_counts) - max(num_nodes - 2, 0);
k_path = zeros(nk_path, 3);
path_dist = zeros(nk_path, 1);
node_locs = zeros(1, num_nodes);
current_dist = 0;
cursor = 1;

for i = 1:(num_nodes - 1)
    k_start = path_nodes(i, :);
    k_stop = path_nodes(i + 1, :);
    n_points = segment_counts(i);
    segment_length = norm(k_stop - k_start);
    t = linspace(0, 1, n_points)';
    if i > 1
        t = t(2:end);
    end
    indices = cursor:(cursor + numel(t) - 1);
    k_path(indices, :) = k_start + t .* (k_stop - k_start);
    path_dist(indices) = current_dist + t * segment_length;
    cursor = indices(end) + 1;
    current_dist = current_dist + segment_length;
    node_locs(i + 1) = current_dist;
end
end

function hr = read_w90_hr(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('compare_w90_bands:CannotOpenHR', ...
        'Cannot open Wannier90 HR file: %s', filename);
end
cleanup_obj = onCleanup(@() fclose(fid));

comment = fgetl(fid); %#ok<NASGU>
num_wann = str2double(strtrim(fgetl(fid)));
nrpts = str2double(strtrim(fgetl(fid)));
if ~isfinite(num_wann) || ~isfinite(nrpts)
    error('compare_w90_bands:BadHRHeader', ...
        'Cannot parse num_wann/nrpts from %s.', filename);
end

degen = zeros(1, 0);
while numel(degen) < nrpts
    line = fgetl(fid);
    if ~ischar(line)
        error('compare_w90_bands:BadHRDegeneracy', ...
            'Unexpected end of file while reading degeneracies.');
    end
    degen = [degen, sscanf(line, '%d').']; %#ok<AGROW>
end
degen = double(degen(1:nrpts));

data = textscan(fid, '%d %d %d %d %d %f %f');
num_rows = numel(data{1});
expected_rows = num_wann * num_wann * nrpts;
if num_rows ~= expected_rows
    error('compare_w90_bands:BadHRBody', ...
        'Expected %d HR matrix rows but found %d in %s.', ...
        expected_rows, num_rows, filename);
end

H_R = complex(zeros(num_wann, num_wann, nrpts));
R_vectors = zeros(nrpts, 3);
block_size = num_wann * num_wann;
for ir = 1:nrpts
    rows = (ir - 1) * block_size + (1:block_size);
    first_row = rows(1);
    R_vectors(ir, :) = [data{1}(first_row), data{2}(first_row), data{3}(first_row)];
    for row = rows
        m = data{4}(row);
        n = data{5}(row);
        H_R(m, n, ir) = complex(data{6}(row), data{7}(row)) / degen(ir);
    end
end

hr = struct();
hr.num_wann = num_wann;
hr.nrpts = nrpts;
hr.degen = degen;
hr.R_vectors = R_vectors;
hr.H_R = H_R;
end

function bands = interpolate_w90_hr(hr, k_path)
nk = size(k_path, 1);
bands = zeros(nk, hr.num_wann);
phase_matrix = exp(1i * 2 * pi * (k_path * hr.R_vectors.'));

for ik = 1:nk
    Hk = complex(zeros(hr.num_wann, hr.num_wann));
    for ir = 1:hr.nrpts
        Hk = Hk + phase_matrix(ik, ir) * hr.H_R(:, :, ir);
    end
    Hk = (Hk + Hk') / 2;
    bands(ik, :) = sort(real(eig(Hk))).';
end
end

function efermi = get_fermi_energy(crystal)
% Crystal objects and struct checkpoints must use the same DFT energy zero.
try
    efermi = crystal.efermi;
catch
    error('compare_w90_bands:MissingFermiEnergy', ...
        'The NSCF checkpoint must provide crystal.efermi in Hartree.');
end
if ~isnumeric(efermi) || ~isscalar(efermi) || ~isreal(efermi) || ~isfinite(efermi)
    error('compare_w90_bands:InvalidFermiEnergy', ...
        'crystal.efermi must be a finite real numeric scalar in Hartree.');
end
end

function [bands, band_indices] = orient_and_select_reference( ...
        bands, reference_size, reference_band_start)
nk_ref = reference_size(1);
nbnd_ref = reference_size(2);
[nrow, ncol] = size(bands);
if nrow == nk_ref
    bands = double(bands);
elseif ncol == nk_ref
    bands = double(bands.');
else
    error('compare_w90_bands:ReferenceKPointMismatch', ...
        ['Reference bands have size [%d, %d], but Wannier90 bands have ', ...
         '%d k-points.'], nrow, ncol, nk_ref);
end

if size(bands, 2) < nbnd_ref
    error('compare_w90_bands:ReferenceBandMismatch', ...
        'Reference has %d bands, but Wannier90 has %d bands.', ...
        size(bands, 2), nbnd_ref);
end

band_start = reference_band_start;
band_stop = band_start + nbnd_ref - 1;
if band_stop > size(bands, 2)
    error('compare_w90_bands:ReferenceBandWindowOutOfRange', ...
        'ReferenceBandStart=%d requires bands %d-%d, but only %d exist.', ...
        band_start, band_start, band_stop, size(bands, 2));
end
band_indices = band_start:band_stop;
bands = bands(:, band_indices);
end

function write_band_csv(filename, path_dist, bands_w90, energy_ref)
[nk, nbnd] = size(bands_w90);
k_index = repelem((1:nk).', nbnd);
band_index = repmat((1:nbnd).', nk, 1);
k_distance = repelem(path_dist(:), nbnd);
w90_eV = reshape(bands_w90.', [], 1);
dft_eV = reshape(energy_ref.', [], 1);
difference_eV = w90_eV - dft_eV;
output_table = table(k_index, k_distance, band_index, ...
    w90_eV, dft_eV, difference_eV);
writetable(output_table, filename);
end
