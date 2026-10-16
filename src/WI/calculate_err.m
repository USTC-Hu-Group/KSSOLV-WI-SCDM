function [err_aver, err_max, err_rms, band_error_detail] = calculate_err( ...
    bands_scdm, energyBands, crystal, varargin)
%CALCULATE_ERR Compare SCDM and DFT band energies on the same k path.
%
% bands_scdm is expected in eV relative to the Fermi level, as returned by
% compute_scdm_only_bands. energyBands may use either [nk x nbnd] or
% [nbnd x nk] orientation.
%
% calculate_err(bands_scdm, energyBands, crystal)
%   Treats energyBands as raw KSSOLV eigenvalues in Hartree and converts
%   them to eV relative to crystal.efermi before comparison.
%
% calculate_err(bands_scdm, energyBands)
%   Treats energyBands as already expressed in eV relative to E_F.
%
% calculate_err(..., 'DiscardTopBands', 4)
%   Excludes the highest 4 interpolated bands from the error metric.
%
% calculate_err(..., 'ReferenceBandStart', 1)
%   When the DFT reference has more bands than the SCDM interpolation,
%   compare against this contiguous DFT band window. The default compares
%   the lowest available DFT bands.
%
% calculate_err(..., 'BandErrorMode', 'entangled', 'Nu', 1, 'Tau', 0.1)
%   Reports Vitale et al. (2020) eta and eta_max. Nu is in eV relative to
%   DFT E_F; Tau is in eV. Neither is read from the SCDM parameters.
%   Isolated eta follows the unnormalized Eq. (18) as printed; err_rms is
%   always the separate ordinary RMS over the selected bands/k-points.

    validate_band_matrix(bands_scdm, 'bands_scdm');
    validate_band_matrix(energyBands, 'energyBands');

    if nargin < 3
        crystal = [];
    elseif is_option_name(crystal)
        varargin = [{crystal}, varargin];
        crystal = [];
    end
    options = parse_error_options(varargin{:});

    if nargin >= 3 && ~isempty(crystal)
        efermi = get_fermi_energy(crystal);
        energyBands = (energyBands - efermi) * 27.211386246;
    end

    [energyBands, reference_band_indices] = orient_and_select_reference_bands( ...
        energyBands, size(bands_scdm), 'energyBands', ...
        options.ReferenceBandStart);

    if ~isreal(bands_scdm) || ~isreal(energyBands)
        error('Band energies must be real-valued.');
    end
    if any(~isfinite(bands_scdm(:))) || any(~isfinite(energyBands(:)))
        error('Band energies contain NaN or Inf values.');
    end

    if options.DiscardTopBands > 0
        if options.DiscardTopBands >= size(bands_scdm, 2)
            error('calculate_err:TooManyDiscardedBands', ...
                ['DiscardTopBands=%d leaves no bands to compare; ', ...
                 'there are only %d bands.'], ...
                options.DiscardTopBands, size(bands_scdm, 2));
        end
        kept = 1:(size(bands_scdm, 2) - options.DiscardTopBands);
        bands_scdm = bands_scdm(:, kept);
        energyBands = energyBands(:, kept);
        reference_band_indices = reference_band_indices(kept);
    end

    difference = bands_scdm - energyBands;
    absolute_difference = abs(difference);
    err_aver = mean(absolute_difference(:));
    err_max = max(absolute_difference(:));
    err_rms = sqrt(mean(difference(:).^2));

    metadata = collect_calculation_metadata(crystal, 'scdm.in');
    [band_eta, band_error_detail] = compute_band_eta( ...
        energyBands, bands_scdm, options);
    band_error_detail.reference_band_indices = reference_band_indices;
    band_error_detail.discarded_top_bands = options.DiscardTopBands;

    fid = fopen('err.txt', 'a');
    if fid == -1
        error('Cannot open err.txt for appending.');
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '%s\n', repmat('=', 1, 72));
    fprintf(fid, 'timestamp %s\n', char(datetime('now', ...
        'Format', 'yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, 'system %s\n', metadata.system_name);
    fprintf(fid, 'working_directory %s\n', pwd);
    fprintf(fid, 'nscf_kgrid %s\n', format_integer_vector(metadata.kgrid));
    fprintf(fid, 'nscf_nkpts %s\n', format_optional_integer(metadata.nkpts));
    fprintf(fid, 'nscf_nbnd %s\n', format_optional_integer(metadata.nbnd));
    fprintf(fid, 'n_wannier %s\n', format_optional_integer(metadata.n_wannier));
    fprintf(fid, 'scdm_filter %s\n', metadata.filter_mode);
    if isfinite(metadata.mu)
        fprintf(fid, 'scdm_mu_eV %.12f\n', metadata.mu);
    end
    if isfinite(metadata.sigma)
        fprintf(fid, 'scdm_sigma_eV %.12f\n', metadata.sigma);
    end
    if isfinite(metadata.scdm_runtime_seconds)
        fprintf(fid, 'scdm_runtime_seconds %.12f\n', ...
            metadata.scdm_runtime_seconds);
        fprintf(fid, 'scdm_runtime_hms %s\n', ...
            format_elapsed_time(metadata.scdm_runtime_seconds));
    else
        fprintf(fid, 'scdm_runtime_seconds unknown\n');
        fprintf(fid, 'scdm_runtime_hms unknown\n');
    end
    fprintf(fid, 'compared_kpoints %d\n', size(bands_scdm, 1));
    fprintf(fid, 'compared_bands %d\n', size(bands_scdm, 2));
    fprintf(fid, 'reference_dft_band_indices %s\n', ...
        format_integer_vector(reference_band_indices));
    fprintf(fid, 'discarded_top_bands %d\n', options.DiscardTopBands);
    fprintf(fid, 'interpolation_method %s\n', options.Method);
    fprintf(fid, 'MAE_eV %.12f\nRMSE_eV %.12f\nMAX_eV %.12f\n', ...
        err_aver, err_rms, err_max);
    fprintf(fid, 'band_error_mode %s\n', band_error_detail.mode);
    fprintf(fid, 'band_error_metric_definition %s\n', band_error_detail.metric_definition);
    fprintf(fid, 'band_error_eta_definition %s\n', band_error_detail.eta_definition);
    fprintf(fid, 'band_error_eta_normalization %s\n', band_error_detail.eta_normalization);
    fprintf(fid, 'band_error_energy_reference DFT_Fermi_level\n');
    fprintf(fid, 'band_error_weight_source %s\n', band_error_detail.weight_source);
    fprintf(fid, 'ETA_%s_eV %.12f\n', band_error_detail.mode, band_eta);
    if isfield(band_error_detail, 'eta_max') && ...
            ~isempty(band_error_detail.eta_max)
        fprintf(fid, 'ETA_MAX_%s_eV %.12f\n', ...
            band_error_detail.mode, band_error_detail.eta_max);
    end
    if strcmp(band_error_detail.mode, 'entangled')
        fprintf(fid, 'band_error_nu_relative_fermi_eV %.17g\n', band_error_detail.nu);
        fprintf(fid, 'band_error_tau_eV %.17g\n', band_error_detail.tau);
    end
    fprintf(fid, 'band_error_eV %s\n', ...
        format_float_vector(band_error_detail.band_error));
    if isfield(band_error_detail, 'band_max_error')
        fprintf(fid, 'band_max_error_eV %s\n', ...
            format_float_vector(band_error_detail.band_max_error));
    end
    fprintf(fid, 'band_weighted_error_eV %s\n', ...
        format_float_vector(band_error_detail.band_weighted_error));
    if isfield(band_error_detail, 'band_weighted_max_error')
        fprintf(fid, 'band_weighted_max_error_eV %s\n', ...
            format_float_vector(band_error_detail.band_weighted_max_error));
    end
    fprintf(fid, 'band_mean_weight %s\n\n', ...
        format_float_vector(band_error_detail.weight));
    clear cleaner;

    fprintf('Band errors: MAE = %.6f eV, RMSE = %.6f eV, MAX = %.6f eV\n', ...
        err_aver, err_rms, err_max);
    fprintf('Band metric eta (%s) = %.6f eV\n', ...
        band_error_detail.mode, band_eta);
    if isfield(band_error_detail, 'eta_max') && ...
            ~isempty(band_error_detail.eta_max)
        fprintf('Band metric eta_max (%s) = %.6f eV\n', ...
            band_error_detail.mode, band_error_detail.eta_max);
    end
    fprintf('Error record appended to %s\n', fullfile(pwd, 'err.txt'));
end

function options = parse_error_options(varargin)
    parser = inputParser;
    parser.PartialMatching = false;
    addParameter(parser, 'DiscardTopBands', 0, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && ...
        isfinite(x) && x >= 0 && x == round(x));
    addParameter(parser, 'Method', 'unknown', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'ReferenceBandStart', 1, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && ...
        isfinite(x) && x >= 1 && x == round(x));
    addParameter(parser, 'BandErrorMode', 'isolated', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'Mu', [], @is_empty_or_scalar_real);
    addParameter(parser, 'Sigma', [], @is_empty_or_scalar_real);
    addParameter(parser, 'Nu', 1.0, ...
        @(x) ~isempty(x) && is_empty_or_scalar_real(x));
    addParameter(parser, 'Tau', 0.1, ...
        @(x) ~isempty(x) && is_empty_or_scalar_real(x) && x > 0);
    addParameter(parser, 'PlotBandError', false, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(parser, 'plot_band_error', [], ...
        @(x) isempty(x) || islogical(x) || isnumeric(x));
    addParameter(parser, 'ErrorPlotPrefix', fullfile('results', 'band_error_metric'), ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'Visible', 'off', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    parse(parser, varargin{:});
    options = parser.Results;
    if ~isempty(options.Mu) || ~isempty(options.Sigma)
        error('calculate_err:SCDMParametersNotEvaluationParameters', ...
            ['Mu/Sigma must not be used for error evaluation. Remove them; ', ...
             'use independent Nu (default 1 eV above E_F) and Tau (0.1 eV).']);
    end
    options.DiscardTopBands = double(options.DiscardTopBands);
    options.Method = char(string(options.Method));
    options.ReferenceBandStart = double(options.ReferenceBandStart);
    options.BandErrorMode = validatestring( ...
        lower(char(string(options.BandErrorMode))), ...
        {'isolated', 'entangled'}, 'calculate_err', 'BandErrorMode');
    if ~isempty(options.plot_band_error)
        options.PlotBandError = logical(options.plot_band_error);
    else
        options.PlotBandError = logical(options.PlotBandError);
    end
    options.ErrorPlotPrefix = char(string(options.ErrorPlotPrefix));
    options.Visible = char(string(options.Visible));
end

function tf = is_option_name(value)
    tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
        any(strcmpi(char(string(value)), ...
        {'DiscardTopBands', 'Method', 'ReferenceBandStart', ...
         'BandErrorMode', 'Mu', 'Sigma', 'Nu', 'Tau', 'PlotBandError', ...
         'plot_band_error', 'ErrorPlotPrefix', 'Visible'}));
end

function tf = is_empty_or_scalar_real(value)
    tf = isempty(value) || ...
        (isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value));
end

function validate_band_matrix(bands, variable_name)
    if ~isnumeric(bands) || ~ismatrix(bands) || isempty(bands)
        error('%s must be a nonempty numeric matrix.', variable_name);
    end
end

function [band_eta, detail] = compute_band_eta( ...
        energyBands, bands_scdm, options)
    metric_args = { ...
        'Nu', options.Nu, 'Tau', options.Tau, ...
        'plot_band_error', options.PlotBandError, ...
        'OutputPrefix', options.ErrorPlotPrefix, ...
        'Visible', options.Visible};

    % calculate_err 内部沿用旧约定 [kpoint x band]；新的模块化
    % calculate_band_error 明确采用 [band x kpoint]。这里统一转置，
    % 避免调用者在外部手动处理矩阵方向。
    [band_eta, detail] = calculate_band_error( ...
        energyBands.', bands_scdm.', options.BandErrorMode, metric_args{:});
end

function [bands, band_indices] = orient_and_select_reference_bands( ...
        bands, reference_size, variable_name, reference_band_start)
    nk_ref = reference_size(1);
    nbnd_ref = reference_size(2);

    [nrow, ncol] = size(bands);
    if nrow == nk_ref
        bands = double(bands);
    elseif ncol == nk_ref
        bands = double(bands.');
    else
        error(['Band k-point dimensions do not match: bands_scdm is ', ...
            '[%d, %d], %s is [%d, %d]. Expected %s to contain %d ', ...
            'k-points in either rows or columns.'], ...
            reference_size(1), reference_size(2), variable_name, ...
            nrow, ncol, variable_name, nk_ref);
    end

    nbnd_available = size(bands, 2);
    if nbnd_available < nbnd_ref
        error(['Band dimensions do not match: bands_scdm has %d bands, ', ...
            'but %s has only %d bands after orientation.'], ...
            nbnd_ref, variable_name, nbnd_available);
    end

    band_start = reference_band_start;
    band_stop = band_start + nbnd_ref - 1;
    if band_stop > nbnd_available
        error('calculate_err:ReferenceBandWindowOutOfRange', ...
            ['ReferenceBandStart=%d with %d SCDM bands requires DFT ', ...
             'bands %d-%d, but %s has only %d bands.'], ...
            band_start, nbnd_ref, band_start, band_stop, ...
            variable_name, nbnd_available);
    end

    band_indices = band_start:band_stop;
    if nbnd_available > nbnd_ref
        warning('calculate_err:ReferenceBandsTrimmed', ...
            ['%s contains %d bands while bands_scdm contains %d bands; ', ...
             'comparing DFT bands %d-%d. Use ReferenceBandStart to ', ...
             'choose a different contiguous DFT band window.'], ...
            variable_name, nbnd_available, nbnd_ref, band_start, band_stop);
    end
    bands = bands(:, band_indices);
end

function efermi = get_fermi_energy(crystal)
    try
        efermi = crystal.efermi;
    catch
        error('The third input must provide a scalar crystal.efermi in Hartree.');
    end

    if ~isnumeric(efermi) || ~isscalar(efermi) || ~isfinite(efermi)
        error('crystal.efermi must be a finite numeric scalar in Hartree.');
    end
end

function metadata = collect_calculation_metadata(crystal, scdm_input)
    metadata = struct( ...
        'system_name', infer_system_name(crystal), ...
        'kgrid', [], ...
        'nkpts', NaN, ...
        'nbnd', NaN, ...
        'n_wannier', NaN, ...
        'filter_mode', 'unknown', ...
        'mu', NaN, ...
        'sigma', NaN, ...
        'scdm_runtime_seconds', NaN);

    if isfile(scdm_input)
        fid = fopen(scdm_input, 'r');
        if fid ~= -1
            cleaner = onCleanup(@() fclose(fid));
            header = sscanf(fgetl(fid), '%g');
            filter_line = strsplit(strtrim(fgetl(fid)));
            clear cleaner;

            if numel(header) >= 3
                metadata.nbnd = header(1);
                metadata.nkpts = header(2);
                metadata.n_wannier = header(3);
            end
            if numel(filter_line) >= 2
                metadata.filter_mode = filter_line{2};
            end
            if numel(filter_line) >= 4
                metadata.mu = str2double(filter_line{3});
                metadata.sigma = str2double(filter_line{4});
            end

            try
                [~, metadata.kgrid] = read_scdm_uniform_kgrid( ...
                    scdm_input, metadata.nkpts);
            catch grid_error
                warning('calculate_err:CannotReadKgrid', ...
                    'Cannot read NSCF KGRID from %s: %s', ...
                    scdm_input, grid_error.message);
            end
        end
    end

    if ~isfinite(metadata.nbnd)
        metadata.nbnd = read_crystal_scalar(crystal, 'nbnd');
    end
    if ~isfinite(metadata.nkpts)
        metadata.nkpts = read_crystal_scalar(crystal, 'nkpts');
    end
    if isempty(metadata.kgrid)
        metadata.kgrid = read_crystal_vector(crystal, 'nkxyz');
    end

    metadata.scdm_runtime_seconds = read_scdm_runtime('scdm_runtime.txt');
end

function system_name = infer_system_name(crystal)
    system_name = read_crystal_text(crystal, 'name');
    if ~isempty(system_name)
        return;
    end

    try
        atoms = crystal.atoms;
        alist = double(crystal.alist(:));
        type_counts = accumarray(alist, 1, [numel(atoms), 1]);
        pieces = {};
        for atom_type = 1:numel(atoms)
            count = type_counts(atom_type);
            if count == 0
                continue;
            end
            symbol = strtrim(char(atoms(atom_type).symbol));
            if count == 1
                pieces{end + 1} = symbol; %#ok<AGROW>
            else
                pieces{end + 1} = sprintf('%s%d', symbol, count); %#ok<AGROW>
            end
        end
        system_name = strjoin(pieces, '');
    catch
        system_name = '';
    end

    if isempty(system_name)
        system_name = 'unknown';
    end
end

function value = read_crystal_scalar(crystal, member)
    value = NaN;
    try
        candidate = crystal.(member);
        if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
            value = double(candidate);
        end
    catch
    end
end

function value = read_crystal_vector(crystal, member)
    value = [];
    try
        candidate = double(crystal.(member));
        if isvector(candidate) && all(isfinite(candidate))
            value = candidate(:).';
        end
    catch
    end
end

function value = read_crystal_text(crystal, member)
    value = '';
    try
        candidate = crystal.(member);
        if ischar(candidate) || (isstring(candidate) && isscalar(candidate))
            value = strtrim(char(candidate));
        end
    catch
    end
end

function text = format_integer_vector(value)
    if isempty(value) || any(~isfinite(value))
        text = 'unknown';
    else
        text = strtrim(sprintf('%d ', round(value)));
    end
end

function text = format_float_vector(value)
    if isempty(value)
        text = 'unknown';
    else
        text = strtrim(sprintf('%.12g ', double(value(:).')));
    end
end

function text = format_optional_integer(value)
    if ~isfinite(value)
        text = 'unknown';
    else
        text = sprintf('%d', round(value));
    end
end

function elapsed_seconds = read_scdm_runtime(filename)
    elapsed_seconds = NaN;
    if ~isfile(filename)
        return;
    end

    try
        text = fileread(filename);
        token = regexp(text, ...
            '(?m)^elapsed_seconds\s+([-+0-9.eE]+)\s*$', ...
            'tokens', 'once');
        if ~isempty(token)
            candidate = str2double(token{1});
            if isfinite(candidate) && candidate >= 0
                elapsed_seconds = candidate;
            end
        end
    catch runtime_error
        warning('calculate_err:CannotReadRuntime', ...
            'Cannot read scdm_m runtime from %s: %s', ...
            filename, runtime_error.message);
    end
end

function text = format_elapsed_time(elapsed_seconds)
    hours = floor(elapsed_seconds / 3600);
    minutes = floor(mod(elapsed_seconds, 3600) / 60);
    seconds = mod(elapsed_seconds, 60);
    text = sprintf('%02d:%02d:%06.3f', hours, minutes, seconds);
end
