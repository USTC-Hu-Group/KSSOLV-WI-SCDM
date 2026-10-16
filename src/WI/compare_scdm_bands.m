function [bands_scdm, H_R, transform_info] = compare_scdm_bands( ...
    kPoints, energyBands, crystal, varargin)
%COMPARE_SCDM_BANDS Plot SCDM-interpolated bands against DFT bands.
%
% kPoints format: {kx, ky, kz, n_points_to_next_node, label; ...}
% The optional first trailing argument can be an AMN filename. Additional
% name-value options are forwarded to compute_scdm_only_bands.

amn_file = 'scdm.amn';
if ~isempty(varargin) && is_text_scalar(varargin{1}) && ...
        ~is_compute_option_name(varargin{1})
    amn_file = char(string(varargin{1}));
    varargin = varargin(2:end);
end
if ~isfile(amn_file)
    error('compare_scdm_bands:MissingAMN', ...
        ['Cannot find AMN file "%s" in the current folder. Pass the same ', ...
         'filename used as the second input of scdm_m.'], amn_file);
end

num_nodes = size(kPoints, 1);
path_nodes = cell2mat(kPoints(:, 1:3));
path_labels = kPoints(:, 5);

segment_counts = cell2mat(kPoints(1:end-1, 4));
if any(segment_counts < 2) || any(segment_counts ~= round(segment_counts))
    error('Each k-path segment must contain an integer number of at least 2 points.');
end

nk_path = sum(segment_counts) - max(num_nodes - 2, 0);
k_path_fine = zeros(nk_path, 3);
path_dist = zeros(nk_path, 1);
node_locs = zeros(1, num_nodes);
current_dist = 0;
node_locs(1) = 0;
cursor = 1;

for i = 1:num_nodes-1
    p_start = path_nodes(i, :);
    p_end = path_nodes(i+1, :);
    n_points = kPoints{i, 4};
    seg_len = norm(p_end - p_start);

    t = linspace(0, 1, n_points)';
    if i > 1
        t = t(2:end);
    end

    seg_kpts = p_start + t .* (p_end - p_start);
    seg_dists = current_dist + t * seg_len;
    segment_indices = cursor:(cursor + numel(t) - 1);

    k_path_fine(segment_indices, :) = seg_kpts;
    path_dist(segment_indices) = seg_dists;
    cursor = segment_indices(end) + 1;

    current_dist = current_dist + seg_len;
    node_locs(i+1) = current_dist;
end

% KSSOLV eband commonly returns [nbnd x nk], while plotting requires
% [nk x nbnd]. Accept either orientation so callers need not transpose.
energyBands = orient_bands_by_kpoints(energyBands, nk_path, 'energyBands');
energybands_ev = (energyBands - crystal.efermi) * 27.211386246;

compute_args = varargin;
if requests_wi_scdm_f(compute_args) && ...
        ~has_compute_option(compute_args, 'TransformEpsilon') && ...
        ~has_compute_option(compute_args, 'TransformParams')
    path_bands_abs_ev = energyBands * 27.211386246;
    finite_path_bands = path_bands_abs_ev(isfinite(path_bands_abs_ev));
    if isempty(finite_path_bands)
        warning('compare_scdm_bands:NoFiniteReferenceBands', ...
            ['Cannot infer Python-compatible TransformEpsilon from ', ...
             'energyBands; falling back to scdm.eig top-band maximum.']);
    else
        python_epsilon = max(finite_path_bands);
        compute_args = [compute_args, {'TransformEpsilon', python_epsilon}];
        fprintf(['WI-SCDM-f Python-compatible epsilon from path bands: ', ...
            '%.12f eV\n'], python_epsilon);
    end
end

[bands_scdm, H_R, transform_info] = compute_scdm_only_bands( ...
    amn_file, 'scdm.eig', 'scdm.in', k_path_fine, crystal.efermi, ...
    compute_args{:});

method_label = format_method_label(transform_info);
method_tag = format_method_tag(transform_info);

fig = figure('Color', 'w', 'Name', ...
    sprintf('Band Structure Comparison (%s)', method_label));
hold on;

plot(path_dist, bands_scdm, 'r--', 'LineWidth', 1.5);
plot(path_dist, energybands_ev, 'k-', 'LineWidth', 1.5);

title(sprintf('%s Interpolated Bands Compared with DFT', method_label));
ylabel('Energy (eV)');
xlim([0, max(path_dist)]);

all_bands = [bands_scdm(:); energybands_ev(:)];
all_bands = all_bands(isfinite(all_bands));
if isempty(all_bands)
    warning('No finite band values found for y-axis limits.');
else
    y_min = min(all_bands);
    y_max = max(all_bands);
    y_pad = max(0.05 * (y_max - y_min), 0.1);
    ylim([y_min - y_pad, y_max + y_pad]);
end

set(gca, 'XTick', node_locs);
set(gca, 'XTickLabel', path_labels);

for loc = node_locs
    xline(loc, 'k:', 'HandleVisibility', 'off');
end

grid on;
box on;
hold off;

if strcmp(method_tag, 'scdm')
    output_file = fullfile(pwd, 'compare_scdm_bands.png');
else
    output_file = fullfile(pwd, sprintf('compare_%s_bands.png', method_tag));
end
exportgraphics(fig, output_file, 'Resolution', 300);
fprintf('Band structure plot saved to %s\n', output_file);

if strcmp(method_tag, 'scdm')
    figure_file = fullfile(pwd, 'compare_scdm_bands.fig');
else
    figure_file = fullfile(pwd, sprintf('compare_%s_bands.fig', method_tag));
end
savefig(fig, figure_file);
fprintf('MATLAB figure saved to %s\n', figure_file);
end

function tf = is_text_scalar(value)
    tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = is_compute_option_name(value)
    value = strtrim(char(string(value)));
    option_names = {'method', 'transformn', 'transforma', ...
        'transformepsilon', 'transformtopband', 'inversesamples', ...
        'transformparams'};
    tf = any(strcmpi(value, option_names));
end

function label = format_method_label(transform_info)
    if isstruct(transform_info) && isfield(transform_info, 'method') && ...
            strcmp(transform_info.method, 'wi-scdm-f')
        label = 'WI-SCDM-f';
    else
        label = 'SCDM-only';
    end
end

function tag = format_method_tag(transform_info)
    if isstruct(transform_info) && isfield(transform_info, 'method') && ...
            strcmp(transform_info.method, 'wi-scdm-f')
        tag = 'wi_scdm_f';
    else
        tag = 'scdm';
    end
end

function bands = orient_bands_by_kpoints(bands, nk, variable_name)
%ORIENT_BANDS_BY_KPOINTS Return a band matrix in [nk x nbnd] orientation.

    if ~isnumeric(bands) || ~ismatrix(bands) || isempty(bands)
        error('%s must be a nonempty numeric matrix.', variable_name);
    end

    [nrow, ncol] = size(bands);
    if nrow == nk
        return;
    end
    if ncol == nk
        bands = bands.';
        return;
    end

    error(['%s dimensions are [%d, %d], but the k path contains %d points. ', ...
        'Expected [nk, nbnd] or [nbnd, nk].'], ...
        variable_name, nrow, ncol, nk);
end

function tf = requests_wi_scdm_f(args)
    tf = false;
    for index = 1:(numel(args) - 1)
        if is_text_scalar(args{index}) && strcmpi(strtrim(char(string(args{index}))), 'Method')
            method = lower(strtrim(char(string(args{index + 1}))));
            tf = any(strcmp(method, {'wi-scdm-f', 'wiscdmf', ...
                'wi_scdm_f', 'scdm-f'}));
            return;
        end
    end
end

function tf = has_compute_option(args, option_name)
    tf = false;
    for index = 1:numel(args)
        if is_text_scalar(args{index}) && ...
                strcmpi(strtrim(char(string(args{index}))), option_name)
            tf = true;
            return;
        end
    end
end
