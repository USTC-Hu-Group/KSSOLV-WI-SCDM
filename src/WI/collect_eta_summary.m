function [summary, pair_summary] = collect_eta_summary(case_root)
%COLLECT_ETA_SUMMARY Collect Vitale 2020 metrics; skip old/other-window eta.
% Entangled: Nu=1 eV relative to DFT E_F, Tau=0.1 eV. Isolated eta
% is the unnormalized printed Eq. (18); RMSE is separate.
% Positive paired delta (MLWF minus SCDM) means MLWF worsens eta.
%
% Usage:
%   collect_eta_summary
%   collect_eta_summary('/path/to/insulator_20_cases')

if nargin < 1 || isempty(case_root)
    case_root = pwd;
end

case_dirs = dir(case_root);
records = cell(0, 19);

for i = 1:numel(case_dirs)
    if ~case_dirs(i).isdir || startsWith(case_dirs(i).name, '.')
        continue;
    end
    err_file = fullfile(case_root, case_dirs(i).name, 'err.txt');
    if ~isfile(err_file)
        continue;
    end
    text = fileread(err_file);
    blocks = regexp(text, '(?m)^={10,}\s*$', 'split');
    for ib = 1:numel(blocks)
        block = strtrim(blocks{ib});
        if ~contains(block, 'interpolation_method')
            continue;
        end
        if ~is_current_vitale_record(block)
            continue;
        end
        method = read_text_field(block, 'interpolation_method');
        eta_mode = read_text_field(block, 'band_error_mode');
        eta_value = read_eta(block);
        eta_max_value = read_eta_max(block, eta_mode);
        records(end + 1, :) = { ...
            case_dirs(i).name, ...
            read_text_field(block, 'system'), ...
            method, ...
            eta_mode, ...
            eta_value, ...
            eta_max_value, ...
            read_numeric_field(block, 'RMSE_eV'), ...
            read_numeric_field(block, 'MAE_eV'), ...
            read_numeric_field(block, 'MAX_eV'), ...
            read_text_field(block, 'n_wannier'), ...
            read_text_field(block, 'nscf_kgrid'), ...
            read_text_field(block, 'nscf_nbnd'), ...
            read_text_field(block, 'timestamp'), ...
            read_text_field(block, 'band_error_metric_definition'), ...
            read_numeric_field(block, 'band_error_nu_relative_fermi_eV'), ...
            read_numeric_field(block, 'band_error_tau_eV'), ...
            read_text_field(block, 'reference_dft_band_indices'), ...
            read_numeric_field(block, 'compared_kpoints'), ...
            read_numeric_field(block, 'compared_bands')}; %#ok<AGROW>
    end
end

summary = cell2table(records, 'VariableNames', { ...
    'case_name', 'system', 'method', 'eta_mode', 'eta_eV', 'eta_max_eV', ...
    'rmse_eV', 'mae_eV', 'max_eV', 'n_wannier', 'nscf_kgrid', ...
    'nscf_nbnd', 'timestamp', 'metric_definition', 'nu_relative_fermi_eV', ...
    'tau_eV', 'reference_dft_band_indices', 'compared_kpoints', 'compared_bands'});

summary_file = fullfile(case_root, 'eta_summary.csv');
writetable(summary, summary_file);
fprintf('Wrote eta summary: %s\n', summary_file);

pair_summary = make_pair_summary(summary);
pair_file = fullfile(case_root, 'eta_pair_summary.csv');
writetable(pair_summary, pair_file);
fprintf('Wrote paired eta summary: %s\n', pair_file);
end

function pair_summary = make_pair_summary(summary)
case_names = unique(summary.case_name, 'stable');
rows = cell(0, 11);
for i = 1:numel(case_names)
    case_name = case_names{i};
    mask = strcmp(summary.case_name, case_name);
    sub = summary(mask, :);
    rows(end + 1, :) = { ...
        case_name, ...
        latest_eta(sub, 'wi-scdm'), ...
        latest_eta_max(sub, 'wi-scdm'), ...
        latest_eta(sub, 'wi-scdm-f'), ...
        latest_eta_max(sub, 'wi-scdm-f'), ...
        latest_eta(sub, 'wlmf'), ...
        latest_eta_max(sub, 'wlmf'), ...
        latest_eta(sub, 'wannier90'), ...
        latest_eta_max(sub, 'wannier90'), ...
        paired_delta(sub, 'wi-scdm'), ...
        paired_delta(sub, 'wi-scdm-f')}; %#ok<AGROW>
end
pair_summary = cell2table(rows, 'VariableNames', { ...
    'case_name', 'wi_scdm_eta_eV', 'wi_scdm_eta_max_eV', ...
    'wi_scdm_f_eta_eV', 'wi_scdm_f_eta_max_eV', ...
    'wlmf_eta_eV', 'wlmf_eta_max_eV', ...
    'wannier90_eta_eV', 'wannier90_eta_max_eV', ...
    'delta_eta_wlmf_minus_wi_scdm_eV', ...
    'delta_eta_wlmf_minus_wi_scdm_f_eV'});
end

function value = latest_eta(table_part, method)
value = NaN;
mask = strcmp(table_part.method, method);
indices = find(mask);
if ~isempty(indices)
    value = table_part.eta_eV(indices(end));
end
end

function value = latest_eta_max(table_part, method)
value = NaN;
mask = strcmp(table_part.method, method);
indices = find(mask);
if ~isempty(indices)
    value = table_part.eta_max_eV(indices(end));
end
end

function value = read_text_field(block, key)
expr = ['(?m)^', regexptranslate('escape', key), '\s+([^\r\n]+)'];
tokens = regexp(block, expr, 'tokens', 'once');
if isempty(tokens)
    value = '';
else
    value = strtrim(tokens{1});
end
end

function value = read_numeric_field(block, key)
text = read_text_field(block, key);
value = str2double(text);
if ~isfinite(value)
    value = NaN;
end
end

function value = read_eta(block)
tokens = regexp(block, '(?m)^ETA_(?:isolated|entangled)_eV\s+([-+0-9.eEdD]+)', ...
    'tokens', 'once');
if isempty(tokens)
    value = NaN;
else
    value = str2double(strrep(tokens{1}, 'D', 'E'));
end
end

function value = read_eta_max(block, eta_mode)
value = NaN;
if ~isempty(eta_mode)
    key = ['ETA_MAX_', eta_mode, '_eV'];
    value = read_numeric_field(block, key);
    if isfinite(value)
        return;
    end
end
tokens = regexp(block, '(?m)^ETA_MAX_(?:isolated|entangled)_eV\s+([-+0-9.eEdD]+)', ...
    'tokens', 'once');
if ~isempty(tokens)
    value = str2double(strrep(tokens{1}, 'D', 'E'));
end
end

function tf = is_current_vitale_record(block)
mode = read_text_field(block, 'band_error_mode');
definition = read_text_field(block, 'band_error_eta_definition');
normalization = read_text_field(block, 'band_error_eta_normalization');
source = read_text_field(block, 'band_error_weight_source');
eta = read_eta(block);
eta_max = read_eta_max(block, mode);
tf = strcmp(read_text_field(block, 'band_error_metric_definition'), 'vitale2020-v1') ...
    && strcmp(read_text_field(block, 'band_error_energy_reference'), 'DFT_Fermi_level') ...
    && isfinite(eta) && eta >= 0 && isfinite(eta_max) && eta_max >= 0;
if strcmp(mode, 'entangled')
    tf = tf && strcmp(definition, 'eq24_fd_geometric_rms') ...
        && strcmp(normalization, 'sum_weights') ...
        && strcmp(source, 'DFT_Wannier_Fermi_Dirac_geometric_mean') ...
        && read_numeric_field(block, 'band_error_nu_relative_fermi_eV') == 1 ...
        && read_numeric_field(block, 'band_error_tau_eV') == 0.1;
else
    tf = tf && strcmp(mode, 'isolated') && strcmp(definition, 'eq18_root_sum_square') ...
        && strcmp(normalization, 'none') && strcmp(source, 'unit_weights');
end
end

function value = paired_delta(sub, method)
value = NaN;
first = find(strcmp(sub.method, 'wlmf') | strcmp(sub.method, 'wannier90') ...
    | strcmp(sub.method, 'mlwf'), 1, 'last');
second = find(strcmp(sub.method, method), 1, 'last');
if isempty(first) || isempty(second)
    return;
end
% Only compare methods on the same reference band set and number of k-points.
if ~strcmp(sub.eta_mode{first}, sub.eta_mode{second}) ...
        || isempty(sub.reference_dft_band_indices{first}) ...
        || ~strcmp(sub.reference_dft_band_indices{first}, sub.reference_dft_band_indices{second}) ...
        || sub.compared_kpoints(first) ~= sub.compared_kpoints(second) ...
        || sub.compared_bands(first) ~= sub.compared_bands(second)
    warning('collect_eta_summary:IncompatiblePair', ...
        'Cannot pair %s: mode, reference bands or k-point counts differ.', sub.case_name{first});
    return;
end
value = sub.eta_eV(first) - sub.eta_eV(second);
end
