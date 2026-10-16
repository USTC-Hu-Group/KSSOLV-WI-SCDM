function [eta, detail] = calculate_band_error(E_DFT, E_Wan, mode, varargin)
%CALCULATE_BAND_ERROR Vitale et al. (2020) DFT/Wannier band distances.
%
% [eta, detail] = calculate_band_error(E_DFT, E_Wan, 'isolated')
% [eta, detail] = calculate_band_error(E_DFT, E_Wan, 'entangled', ...
%     'Nu', 1.0, 'Tau', 0.1)
%
% Input convention:
%   E_DFT and E_Wan must be numeric matrices with size [band x kpoint].
%   Both energy matrices MUST be in eV relative to the SAME DFT Fermi level.
%   Nu is relative to that Fermi level, not the SCDM chemical potential.
%   Isolated input must contain the matched target occupied manifold only.
%   DOI: 10.1038/s41524-020-0312-y, Eqs. (18), (19), (24)-(26).
%   Eq. (18) as printed is an UNNORMALIZED root sum of squares. The ordinary
%   RMS is returned separately as detail.rmse; eta must not be labelled RMS
%   for isolated bands. Entangled eta is normalized by the weight sum.
%
% Output:
%   eta                         Eq. (18) or Eq. (24), in eV
%   detail.eta_max              max band-distance metric
%   detail.band_error           per-band ordinary RMS error
%   detail.band_max_error       per-band ordinary max absolute error
%   detail.band_weighted_error  per-band weighted RMS error
%   detail.band_weighted_max_error
%                               per-band weighted max absolute error
%   detail.weight               per-band mean weight

    if nargin < 3
        error('calculate_band_error:NotEnoughInputs', ...
            'E_DFT, E_Wan, and mode are required.');
    end

    parser = inputParser;
    parser.FunctionName = mfilename;
    parser.PartialMatching = false;
    addParameter(parser, 'Nu', 1.0, @is_scalar_real);
    addParameter(parser, 'Tau', 0.1, @(x) is_scalar_real(x) && x > 0);
    % Accept empty legacy placeholders, but never reinterpret SCDM inputs.
    addParameter(parser, 'mu', [], @is_empty_or_scalar_real);
    addParameter(parser, 'sigma', [], @is_empty_or_scalar_real);
    addParameter(parser, 'plot_band_error', false, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(parser, 'PlotBandError', [], ...
        @(x) isempty(x) || islogical(x) || isnumeric(x));
    addParameter(parser, 'OutputPrefix', '', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'Visible', 'off', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    parse(parser, varargin{:});
    if ~isempty(parser.Results.mu) || ~isempty(parser.Results.sigma)
        error('calculate_band_error:SCDMParametersNotEvaluationParameters', ...
            ['Mu/Sigma are SCDM construction parameters, not band-metric ', ...
             'parameters. Remove them; use Nu (default 1 eV above E_F) ', ...
             'and Tau (default 0.1 eV) for the Vitale metric.']);
    end

    mode = validatestring(lower(char(string(mode))), ...
        {'isolated', 'entangled'}, mfilename, 'mode');
    E_DFT = validate_energy_matrix(E_DFT, 'E_DFT');
    E_Wan = validate_energy_matrix(E_Wan, 'E_Wan');

    if ~isequal(size(E_DFT), size(E_Wan))
        error('calculate_band_error:SizeMismatch', ...
            ['E_DFT and E_Wan must have the same size [band x kpoint]. ', ...
             'Got E_DFT=[%d,%d], E_Wan=[%d,%d].'], ...
            size(E_DFT, 1), size(E_DFT, 2), size(E_Wan, 1), size(E_Wan, 2));
    end

    [nbnd, nkpts] = size(E_DFT);
    difference = E_DFT - E_Wan;
    absolute_difference = abs(difference);
    band_error = sqrt(mean(difference.^2, 2));
    band_max_error = max(absolute_difference, [], 2);

    nu = double(parser.Results.Nu);
    tau = double(parser.Results.Tau);

    switch mode
        case 'isolated'
            % 孤立能带具有清晰的一一对应关系，因此每个 band、每个 k 点
            % 都应同等参与误差统计，不引入任何能量或费米权重。
            weight_matrix = ones(nbnd, nkpts);
            weighted_absolute_difference = absolute_difference;
            band_mean_weight = ones(nbnd, 1);
            band_weighted_error = band_error;
            band_weighted_max_error = band_max_error;
            eta = norm(difference(:)); % Eq. (18), exactly as printed.
            % Vitale et al. Eq. (19): eta_max selects the single
            % worst-interpolated band/k-point for an isolated manifold.
            eta_max = max(absolute_difference(:));
            nu = [];
            tau = [];
            eta_definition = 'eq18_root_sum_square';
            eta_normalization = 'none';
            weight_source = 'unit_weights';

        case 'entangled'
            % Eq. (25): geometric mean of TWO Fermi-Dirac occupations.
            % Evaluate in log space to avoid exp overflow and premature
            % underflow of f_DFT*f_Wan. Never reuse the SCDM erfc filter.
            log_f_dft = log_fermi_dirac((E_DFT - nu) ./ tau);
            log_f_wan = log_fermi_dirac((E_Wan - nu) ./ tau);
            weight_matrix = exp(0.5 * log_f_dft + 0.5 * log_f_wan);
            total_weight = sum(weight_matrix(:));
            if total_weight <= 0
                error('calculate_band_error:ZeroWeight', ...
                    ['The Fermi-Dirac weight sum is numerically zero. ', ...
                     'Check Nu/Tau and the common DFT Fermi-energy zero.']);
            end

            eta = norm(sqrt(weight_matrix(:) / total_weight) .* difference(:));
            % Eq. (26): weight*abs(delta), NOT sqrt(weight)*abs(delta).
            weighted_absolute_difference = weight_matrix .* absolute_difference;
            eta_max = max(weighted_absolute_difference(:));
            band_weight_sum = sum(weight_matrix, 2);
            band_mean_weight = band_weight_sum / nkpts;
            band_weighted_error = NaN(nbnd, 1);
            valid_band = band_weight_sum > 0;
            band_weighted_error(valid_band) = sqrt( ...
                sum(weight_matrix(valid_band, :) .* ...
                difference(valid_band, :).^2, 2) ./ band_weight_sum(valid_band));
            if any(~valid_band)
                warning('calculate_band_error:ZeroBandWeight', ...
                    'Some bands have zero total weight; their weighted band error is NaN.');
            end
            band_weighted_max_error = max(weighted_absolute_difference, [], 2);
            eta_definition = 'eq24_fd_geometric_rms';
            eta_normalization = 'sum_weights';
            weight_source = 'DFT_Wannier_Fermi_Dirac_geometric_mean';
    end

    detail = struct();
    detail.mode = mode;
    detail.metric_definition = 'vitale2020-v1';
    detail.eta_definition = eta_definition;
    detail.eta_normalization = eta_normalization;
    detail.weight_source = weight_source;
    detail.energy_reference = 'DFT_Fermi_level';
    detail.mae = mean(absolute_difference(:));
    detail.rmse = norm(difference(:)) / sqrt(numel(difference));
    detail.max_absolute_error = max(absolute_difference(:));
    detail.eta = eta;
    detail.eta_max = eta_max;
    detail.band_error = band_error;
    detail.band_max_error = band_max_error;
    detail.band_weighted_error = band_weighted_error;
    detail.band_weighted_max_error = band_weighted_max_error;
    detail.weight = band_mean_weight;
    detail.weight_matrix = weight_matrix;
    detail.weighted_absolute_difference = weighted_absolute_difference;
    detail.nu = nu;
    detail.tau = tau;
    detail.nbnd = nbnd;
    detail.nkpts = nkpts;
    detail.difference = difference;
    detail.plot_files = struct();

    plot_band_error = logical(parser.Results.plot_band_error);
    if ~isempty(parser.Results.PlotBandError)
        plot_band_error = logical(parser.Results.PlotBandError);
    end
    if plot_band_error
        detail.plot_files = plot_error_diagnostics( ...
            E_DFT, mode, detail, char(parser.Results.OutputPrefix), ...
            char(parser.Results.Visible));
    end
end

function value = log_fermi_dirac(x)
    value = -max(x, 0) - log1p(exp(-abs(x)));
end

function tf = is_scalar_real(value)
    tf = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end

function matrix = validate_energy_matrix(matrix, variable_name)
    if ~isnumeric(matrix) || ~ismatrix(matrix) || isempty(matrix)
        error('calculate_band_error:BadInput', ...
            '%s must be a nonempty numeric matrix [band x kpoint].', ...
            variable_name);
    end
    if ~isreal(matrix)
        error('calculate_band_error:ComplexEnergy', ...
            '%s must contain real-valued band energies.', variable_name);
    end
    matrix = double(matrix);
    if any(~isfinite(matrix(:)))
        error('calculate_band_error:NonFiniteEnergy', ...
            '%s contains NaN or Inf values.', variable_name);
    end
end

function tf = is_empty_or_scalar_real(value)
    tf = isempty(value) || ...
        (isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value));
end

function plot_files = plot_error_diagnostics(E_DFT, mode, detail, output_prefix, visible)
    if isempty(output_prefix)
        output_prefix = fullfile(pwd, ['band_error_', mode]);
    end
    output_dir = fileparts(output_prefix);
    if ~isempty(output_dir) && ~isfolder(output_dir)
        mkdir(output_dir);
    end

    plot_files = struct();
    band_index = (1:detail.nbnd).';
    fig = figure('Color', 'w', 'Visible', visible, ...
        'Name', ['Band Error - ', mode]);
    if strcmp(mode, 'entangled')
        bar(band_index, detail.band_weighted_error, 0.75);
        ylabel('Weighted RMS error (eV)');
        title('Entangled Band Error by Band');
    else
        bar(band_index, detail.band_error, 0.75);
        ylabel('RMS error (eV)');
        title('Isolated Band Error by Band');
    end
    xlabel('Band index');
    grid on;
    box on;
    plot_files.band_error_png = [output_prefix, '_band_error.png'];
    plot_files.band_error_fig = [output_prefix, '_band_error.fig'];
    save_figure_pair(fig, plot_files.band_error_png, plot_files.band_error_fig);
    close(fig);

    if strcmp(mode, 'entangled')
        fig = figure('Color', 'w', 'Visible', visible, ...
            'Name', 'Band Energy vs Vitale Evaluation Weight');
        band_ids = repmat((1:detail.nbnd).', 1, detail.nkpts);
        scatter(E_DFT(:), detail.weight_matrix(:), 18, band_ids(:), 'filled');
        xlabel('DFT energy - E_F (eV)');
        ylabel('Geometric-mean Fermi-Dirac weight');
        title(sprintf('Vitale weight: nu = %.3g eV, tau = %.3g eV', ...
            detail.nu, detail.tau));
        cb = colorbar;
        cb.Label.String = 'Band index';
        ylim([-0.02, 1.02]);
        grid on;
        box on;
        plot_files.energy_weight_png = [output_prefix, '_energy_weight.png'];
        plot_files.energy_weight_fig = [output_prefix, '_energy_weight.fig'];
        save_figure_pair(fig, plot_files.energy_weight_png, ...
            plot_files.energy_weight_fig);
        close(fig);
    end
end

function save_figure_pair(fig, png_file, fig_file)
    try
        exportgraphics(fig, png_file, 'Resolution', 300);
    catch
        print(fig, png_file, '-dpng', '-r300');
    end
    savefig(fig, fig_file);
end
