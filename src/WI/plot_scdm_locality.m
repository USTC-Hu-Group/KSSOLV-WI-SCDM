function [output_file, locality] = plot_scdm_locality(varargin)
%PLOT_SCDM_LOCALITY Quantify SCDM orbital locality from Cube files.
%
% [output_file, locality] = plot_scdm_locality()
% [output_file, locality] = plot_scdm_locality(logspace(-1, -4, 41))
% [output_file, locality] = plot_scdm_locality('Epsilons', logspace(-1, -4, 41))
%
% The function reads all *.cube files in the current folder by default.
% For each relative truncation threshold epsilon, it counts the fraction of
% entries that remain nonzero after setting entries with
% |psi(r)| < epsilon * max(|psi(r)|) to zero. The fraction is averaged over
% all localized functions, converted to nz%, and plotted as (nz%)^(1/3)
% versus log10(epsilon). Exponential localization appears as an
% approximately straight line.

    if ~isempty(varargin) && isnumeric(varargin{1})
        varargin = [{'Epsilons'}, varargin];
    end

    parser = inputParser;
    addParameter(parser, 'CubePattern', '*.cube', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'Epsilons', logspace(-1, -3, 31), ...
        @(x) isnumeric(x) && isvector(x) && all(isfinite(x(:))) && ...
        all(x(:) > 0) && all(x(:) <= 1));
    addParameter(parser, 'OutputPrefix', 'scdm_locality', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'Resolution', 450, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'MaxFiles', Inf, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'LinearR2Target', 0.987, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    addParameter(parser, 'MinFitPoints', 6, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(parser, 'TargetVisualSlope', -1, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x ~= 0);
    addParameter(parser, 'ShowDiscarded', false, ...
        @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(parser, 'Visible', 'off', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    parse(parser, varargin{:});
    options = parser.Results;

    cube_files = dir(char(options.CubePattern));
    if isempty(cube_files)
        error('plot_scdm_locality:NoCubeFiles', ...
            'No Cube files matching "%s" were found in %s.', ...
            char(options.CubePattern), pwd);
    end
    [~, order] = sort({cube_files.name});
    cube_files = cube_files(order);
    if isfinite(options.MaxFiles)
        cube_files = cube_files(1:min(numel(cube_files), options.MaxFiles));
    end

    epsilons = options.Epsilons(:);
    epsilons = sort(epsilons, 'descend');
    epsilons = epsilons([true; diff(epsilons) ~= 0]);
    n_eps = numel(epsilons);
    n_orb = numel(cube_files);

    nz_fraction_per_orbital = NaN(n_eps, n_orb);
    grid_points = NaN(n_orb, 1);
    max_abs_values = NaN(n_orb, 1);
    filenames = strings(n_orb, 1);

    for iorb = 1:n_orb
        filenames(iorb) = string(cube_files(iorb).name);
        [values, grid_size] = read_cube_values(cube_files(iorb).name);
        amplitudes = abs(values(:));
        max_abs = max(amplitudes);
        max_abs_values(iorb) = max_abs;
        grid_points(iorb) = prod(grid_size);
        if max_abs <= 0 || ~isfinite(max_abs)
            warning('plot_scdm_locality:ZeroCube', ...
                'Cube file %s has zero or nonfinite maximum amplitude.', ...
                cube_files(iorb).name);
            continue;
        end

        for ieps = 1:n_eps
            cutoff = epsilons(ieps) * max_abs;
            nz_fraction_per_orbital(ieps, iorb) = ...
                nnz(amplitudes >= cutoff) / numel(amplitudes);
        end
    end

    nz_fraction = mean(nz_fraction_per_orbital, 2, 'omitnan');
    nz_percent = 100 * nz_fraction;
    diameter_metric = nz_percent.^(1/3);
    log10_epsilon = log10(epsilons);

    linear_mask = select_linear_region(log10_epsilon, diameter_metric, options);
    discarded_mask = ~linear_mask & isfinite(log10_epsilon) & ...
        isfinite(diameter_metric);
    [fit_y, fit_slope, fit_intercept, fit_r2, fit_C] = fit_linear_logeps( ...
        log10_epsilon(linear_mask), diameter_metric(linear_mask), ...
        log10_epsilon(linear_mask));

    output_file = make_output_filename(options.OutputPrefix, n_orb, epsilons);
    fig = figure('Color', 'w', 'Visible', char(options.Visible), ...
        'Units', 'pixels', 'Position', [100, 100, 900, 650]);
    ax = axes(fig);
    hold(ax, 'on');

    plot_color = [0.02, 0.09, 0.36];
    legend_handles = gobjects(0);
    legend_labels = {};

    if options.ShowDiscarded && any(discarded_mask)
        [x_discarded, discarded_order] = sort(log10_epsilon(discarded_mask), ...
            'ascend');
        y_discarded = diameter_metric(discarded_mask);
        y_discarded = y_discarded(discarded_order);
        h_discarded = plot(ax, x_discarded, y_discarded, 'o', ...
            'Color', [0.70, 0.70, 0.70], ...
            'MarkerFaceColor', [0.86, 0.86, 0.86], ...
            'MarkerEdgeColor', [0.70, 0.70, 0.70], ...
            'MarkerSize', 5);
        legend_handles(end + 1) = h_discarded;
        legend_labels{end + 1} = 'Discarded tail';
    end

    [x_selected, selected_order] = sort(log10_epsilon(linear_mask), 'ascend');
    y_selected = diameter_metric(linear_mask);
    y_selected = y_selected(selected_order);
    h_data = plot(ax, x_selected, y_selected, '-o', ...
        'Color', plot_color, ...
        'MarkerFaceColor', [0.96, 0.74, 0.12], ...
        'MarkerEdgeColor', plot_color, ...
        'LineWidth', 1.8, ...
        'MarkerSize', 6);
    legend_handles(end + 1) = h_data;
    legend_labels{end + 1} = 'SCDM Cube data';

    if ~isempty(fit_y)
        fit_y = fit_y(selected_order);
        h_fit = plot(ax, x_selected, fit_y, '--', ...
            'Color', [0.75, 0.10, 0.08], ...
            'LineWidth', 1.5);
        legend_handles(end + 1) = h_fit;
        legend_labels{end + 1} = sprintf( ...
            'linear fit: slope = %.3g, R^2 = %.4f', fit_slope, fit_r2);
    end
    legend(ax, legend_handles, legend_labels, ...
        'Location', 'southwest', 'Box', 'off');

    xlabel(ax, '$\log_{10}(\epsilon)$', 'Interpreter', 'latex');
    ylabel(ax, '$(nz\%)^{1/3}$', 'Interpreter', 'latex');
    title(ax, sprintf('SCDM locality from %d Cube orbitals', n_orb), ...
        'FontName', 'Arial', 'FontWeight', 'bold');
    grid(ax, 'on');
    ax.Box = 'on';
    ax.FontName = 'Arial';
    ax.FontSize = 12;
    ax.LineWidth = 1.0;
    set_axis_limits(ax, x_selected, y_selected);
    adjust_plot_aspect_for_slope(ax, fit_slope, options.TargetVisualSlope);

    exportgraphics(fig, output_file, 'Resolution', options.Resolution);

    locality = struct();
    locality.files = filenames;
    locality.epsilon = epsilons;
    locality.log10_epsilon = log10_epsilon;
    locality.nz_fraction = nz_fraction;
    locality.nz_percent = nz_percent;
    locality.diameter_metric = diameter_metric;
    locality.linear_mask = linear_mask;
    locality.discarded_mask = discarded_mask;
    locality.linear_epsilon = epsilons(linear_mask);
    locality.linear_log10_epsilon = log10_epsilon(linear_mask);
    locality.linear_diameter_metric = diameter_metric(linear_mask);
    locality.nz_fraction_per_orbital = nz_fraction_per_orbital;
    locality.grid_points = grid_points;
    locality.max_abs = max_abs_values;
    locality.fit_C = fit_C;
    locality.fit_intercept = fit_intercept;
    locality.fit_slope_log10 = fit_slope;
    locality.fit_r2 = fit_r2;
    locality.output_file = fullfile(pwd, output_file);

    fprintf('SCDM locality image written to %s\n', locality.output_file);
    fprintf('Linear region: %d/%d points kept; discarded %d low-threshold tail points.\n', ...
        nnz(linear_mask), numel(linear_mask), nnz(discarded_mask));
    fprintf('Fit model: (nz%%)^(1/3) ~= %.6g * log10(epsilon) + %.6g, R^2 = %.6f\n', ...
        fit_slope, fit_intercept, fit_r2);
    fprintf('Equivalent exp decay constant C ~= %.6g\n', fit_C);

end

function [values, grid_size] = read_cube_values(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('plot_scdm_locality:CannotOpenCube', ...
            'Cannot open Cube file %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    fgetl(fid);
    fgetl(fid);
    atom_line = sscanf(fgetl(fid), '%f');
    if numel(atom_line) < 4
        error('plot_scdm_locality:BadCubeHeader', ...
            'Cube file %s has an invalid atom/origin line.', filename);
    end
    natoms = abs(round(atom_line(1)));

    grid_size = zeros(1, 3);
    for axis_index = 1:3
        line_values = sscanf(fgetl(fid), '%f');
        if numel(line_values) < 4
            error('plot_scdm_locality:BadCubeGrid', ...
                'Cube file %s has an invalid grid line.', filename);
        end
        grid_size(axis_index) = abs(round(line_values(1)));
    end

    for atom_index = 1:natoms
        atom_values = sscanf(fgetl(fid), '%f');
        if numel(atom_values) < 5
            error('plot_scdm_locality:BadCubeAtoms', ...
                'Cube file %s has an invalid atom line.', filename);
        end
    end

    expected_values = prod(grid_size);
    raw_values = fscanf(fid, '%f');
    if numel(raw_values) < expected_values
        error('plot_scdm_locality:BadCubeDataLength', ...
            'Cube file %s contains %d values; expected at least %d.', ...
            filename, numel(raw_values), expected_values);
    elseif numel(raw_values) > expected_values
        warning('plot_scdm_locality:ExtraCubeData', ...
            'Cube file %s contains extra values; using the first data block.', ...
            filename);
        raw_values = raw_values(1:expected_values);
    end
    values = raw_values(:);
end

function linear_mask = select_linear_region(log10_epsilon, diameter_metric, options)
    valid_mask = isfinite(log10_epsilon) & isfinite(diameter_metric) & ...
        diameter_metric > 0;
    valid_indices = find(valid_mask);
    linear_mask = false(size(log10_epsilon));
    if isempty(valid_indices)
        return;
    end

    min_points = min(numel(valid_indices), max(2, round(options.MinFitPoints)));
    best_end = min_points;
    best_score = -Inf;
    found_target = false;

    for endpoint = numel(valid_indices):-1:min_points
        candidate_indices = valid_indices(1:endpoint);
        [slope, ~, r2] = linear_fit_stats( ...
            log10_epsilon(candidate_indices), diameter_metric(candidate_indices));
        if slope >= 0 || ~isfinite(r2)
            continue;
        end
        score = r2 + 0.001 * endpoint;
        if score > best_score
            best_score = score;
            best_end = endpoint;
        end
        if r2 >= options.LinearR2Target
            best_end = endpoint;
            found_target = true;
            break;
        end
    end

    if ~found_target && best_score == -Inf
        best_end = min_points;
    end
    linear_mask(valid_indices(1:best_end)) = true;
end

function [fit_y, fit_slope, fit_intercept, fit_r2, fit_C] = ...
        fit_linear_logeps(log10_epsilon, diameter_metric, query_log10_epsilon)
    fit_y = [];
    fit_slope = NaN;
    fit_intercept = NaN;
    fit_r2 = NaN;
    fit_C = NaN;
    if numel(log10_epsilon) < 2
        return;
    end

    [fit_slope, fit_intercept, fit_r2] = linear_fit_stats( ...
        log10_epsilon, diameter_metric);
    fit_y = fit_slope * query_log10_epsilon(:) + fit_intercept;
    if fit_slope < 0
        fit_C = -log(10) / fit_slope;
    end
end

function [slope, intercept, r2] = linear_fit_stats(x, y)
    coefficients = polyfit(x(:), y(:), 1);
    slope = coefficients(1);
    intercept = coefficients(2);
    predicted = polyval(coefficients, x(:));
    observed = y(:);
    residual_sum = sum((observed - predicted).^2);
    total_sum = sum((observed - mean(observed)).^2);
    r2 = NaN;
    if total_sum > 0
        r2 = 1 - residual_sum / total_sum;
    end
end

function set_axis_limits(ax, x_data, y_data)
    if isempty(x_data) || isempty(y_data)
        axis(ax, 'tight');
        return;
    end

    x_min = min(x_data);
    x_max = max(x_data);
    y_min = min(y_data);
    y_max = max(y_data);
    x_pad = 0.06 * max(eps, x_max - x_min);
    y_pad = 0.08 * max(eps, y_max - y_min);
    xlim(ax, [x_min - x_pad, x_max + x_pad]);
    ylim(ax, [y_min - y_pad, y_max + y_pad]);
end

function adjust_plot_aspect_for_slope(ax, fit_slope, target_visual_slope)
    if ~isfinite(fit_slope) || fit_slope == 0
        return;
    end

    x_limits = xlim(ax);
    y_limits = ylim(ax);
    range_x = diff(x_limits);
    range_y = diff(y_limits);
    if range_x <= 0 || range_y <= 0
        return;
    end

    normalized_slope = fit_slope * range_x / range_y;
    if ~isfinite(normalized_slope) || normalized_slope == 0
        return;
    end

    height_over_width = abs(target_visual_slope / normalized_slope);
    height_over_width = min(max(height_over_width, 0.45), 2.20);
    pbaspect(ax, [1, height_over_width, 1]);
end

function output_file = make_output_filename(prefix, n_orb, epsilons)
    eps_min = min(epsilons);
    eps_max = max(epsilons);
    output_file = sprintf('%s_%dorb_eps_%s_to_%s.png', char(prefix), n_orb, ...
        format_epsilon_for_file(eps_min), format_epsilon_for_file(eps_max));
end

function text = format_epsilon_for_file(value)
    text = sprintf('%.0e', value);
    text = strrep(text, '+', '');
    text = strrep(text, '-', 'm');
end
