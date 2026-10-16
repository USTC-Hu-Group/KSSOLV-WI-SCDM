function tf = has_vitale_band_error(err_file, method, mode, nu, tau)
%HAS_VITALE_BAND_ERROR Validate the latest exact-method metric record.
% Records from the erfc evaluator, other energy windows, or partial writes
% must be recomputed. In particular, wi-scdm must not match wi-scdm-f.
tf = false;
if ~isfile(err_file)
    return;
end
blocks = regexp(fileread(err_file), '(?m)^={10,}[^\S\r\n]*$', 'split');
for i = numel(blocks):-1:1
    block = blocks{i};
    if ~strcmp(read_field(block, 'interpolation_method'), method)
        continue;
    end
    tf = strcmp(read_field(block, 'band_error_metric_definition'), 'vitale2020-v1') ...
        && strcmp(read_field(block, 'band_error_mode'), mode) ...
        && strcmp(read_field(block, 'band_error_energy_reference'), 'DFT_Fermi_level') ...
        && is_valid_metric(block, ['ETA_', mode, '_eV']) ...
        && is_valid_metric(block, ['ETA_MAX_', mode, '_eV']);
    if strcmp(mode, 'entangled')
        tf = tf && strcmp(read_field(block, 'band_error_eta_definition'), ...
            'eq24_fd_geometric_rms') ...
            && strcmp(read_field(block, 'band_error_eta_normalization'), 'sum_weights') ...
            && strcmp(read_field(block, 'band_error_weight_source'), ...
            'DFT_Wannier_Fermi_Dirac_geometric_mean') ...
            && str2double(read_field(block, 'band_error_nu_relative_fermi_eV')) == nu ...
            && str2double(read_field(block, 'band_error_tau_eV')) == tau;
    else
        tf = tf && strcmp(read_field(block, 'band_error_eta_definition'), ...
            'eq18_root_sum_square') ...
            && strcmp(read_field(block, 'band_error_eta_normalization'), 'none') ...
            && strcmp(read_field(block, 'band_error_weight_source'), 'unit_weights');
    end
    return;
end
end

function value = read_field(block, key)
token = regexp(block, ['(?m)^', regexptranslate('escape', key), ...
    '[^\S\r\n]+([^\r\n]+)'], 'tokens', 'once');
value = '';
if ~isempty(token)
    value = strtrim(token{1});
end
end

function tf = is_valid_metric(block, key)
value = str2double(read_field(block, key));
tf = isfinite(value) && value >= 0;
end
