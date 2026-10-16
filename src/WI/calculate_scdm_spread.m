function [spread_table, summary] = calculate_scdm_spread(varargin)
%CALCULATE_SCDM_SPREAD Compute real-space spreads of SCDM Wannier orbitals.
%
% Usage from Cube files in the current folder:
%   [spread_table, summary] = calculate_scdm_spread()
%   [spread_table, summary] = calculate_scdm_spread('CubePattern', 'si_*.cube')
%
% Usage from in-memory SCDM orbitals:
%   [spread_table, summary] = calculate_scdm_spread(wan_r, grid_size, structure)
%
% The spread of orbital n is evaluated as
%   Omega_n = int |psi_n(r)|^2 |r - r_n|^2 dr / int |psi_n(r)|^2 dr.
% For periodic systems the default center and displacement are computed in
% fractional coordinates using circular means and minimum-image wrapping.

    bohr_to_ang = 0.529177210903;

    [source, options] = parse_inputs(varargin{:});
    orbitals = load_orbitals(source);
    n_orb = numel(orbitals);
    if n_orb == 0
        error('calculate_scdm_spread:NoOrbitals', ...
            'No SCDM orbitals were found.');
    end

    rows = repmat(empty_spread_row(), n_orb, 1);
    for iorb = 1:n_orb
        rows(iorb) = compute_one_spread(orbitals(iorb), options.Method);
    end

    orbital_index = (1:n_orb).';
    filename = strings(n_orb, 1);
    norm_value = zeros(n_orb, 1);
    spread_bohr2 = zeros(n_orb, 1);
    spread_ang2 = zeros(n_orb, 1);
    radius_bohr = zeros(n_orb, 1);
    radius_ang = zeros(n_orb, 1);
    center_frac = zeros(n_orb, 3);
    center_bohr = zeros(n_orb, 3);
    center_ang = zeros(n_orb, 3);
    max_abs = zeros(n_orb, 1);
    grid_points = zeros(n_orb, 1);

    for iorb = 1:n_orb
        filename(iorb) = string(rows(iorb).filename);
        norm_value(iorb) = rows(iorb).norm_value;
        spread_bohr2(iorb) = rows(iorb).spread_bohr2;
        spread_ang2(iorb) = rows(iorb).spread_bohr2 * bohr_to_ang^2;
        radius_bohr(iorb) = sqrt(max(rows(iorb).spread_bohr2, 0));
        radius_ang(iorb) = radius_bohr(iorb) * bohr_to_ang;
        center_frac(iorb, :) = rows(iorb).center_frac;
        center_bohr(iorb, :) = rows(iorb).center_bohr;
        center_ang(iorb, :) = rows(iorb).center_bohr * bohr_to_ang;
        max_abs(iorb) = rows(iorb).max_abs;
        grid_points(iorb) = rows(iorb).grid_points;
    end

    spread_table = table(orbital_index, filename, norm_value, ...
        spread_bohr2, spread_ang2, radius_bohr, radius_ang, ...
        center_frac(:, 1), center_frac(:, 2), center_frac(:, 3), ...
        center_bohr(:, 1), center_bohr(:, 2), center_bohr(:, 3), ...
        center_ang(:, 1), center_ang(:, 2), center_ang(:, 3), ...
        max_abs, grid_points, ...
        'VariableNames', {'orbital', 'file', 'norm', ...
        'spread_bohr2', 'spread_ang2', 'rms_radius_bohr', ...
        'rms_radius_ang', 'center_f1', 'center_f2', 'center_f3', ...
        'center_x_bohr', 'center_y_bohr', 'center_z_bohr', ...
        'center_x_ang', 'center_y_ang', 'center_z_ang', ...
        'max_abs_value', 'grid_points'});

    summary = struct();
    summary.n_orbitals = n_orb;
    summary.method = options.Method;
    summary.total_spread_bohr2 = sum(spread_bohr2, 'omitnan');
    summary.total_spread_ang2 = sum(spread_ang2, 'omitnan');
    summary.mean_spread_bohr2 = mean(spread_bohr2, 'omitnan');
    summary.mean_spread_ang2 = mean(spread_ang2, 'omitnan');
    summary.min_spread_ang2 = min(spread_ang2);
    summary.max_spread_ang2 = max(spread_ang2);
    summary.source = source.kind;

    if options.WriteFiles
        [summary.text_file, summary.csv_file] = write_spread_outputs( ...
            spread_table, summary, options);
    else
        summary.text_file = "";
        summary.csv_file = "";
    end

    print_summary(summary);
end

function [source, options] = parse_inputs(varargin)
    source = struct('kind', 'cube', 'cube_pattern', '*.cube', ...
        'wan_r', [], 'grid_size', [], 'structure', []);

    if ~isempty(varargin) && isnumeric(varargin{1})
        if numel(varargin) < 3
            error('calculate_scdm_spread:MissingDirectInputs', ...
                ['Direct mode requires calculate_scdm_spread(wan_r, ', ...
                 'grid_size, structure, ...).']);
        end
        source.kind = 'direct';
        source.wan_r = varargin{1};
        source.grid_size = varargin{2};
        source.structure = varargin{3};
        varargin = varargin(4:end);
    end

    parser = inputParser;
    addParameter(parser, 'CubePattern', '*.cube', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'Method', 'periodic', ...
        @(x) any(strcmpi(char(string(x)), {'periodic', 'direct'})));
    addParameter(parser, 'OutputPrefix', 'scdm_spread', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'WriteFiles', true, ...
        @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(parser, 'Overwrite', false, ...
        @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    parse(parser, varargin{:});

    options = parser.Results;
    options.CubePattern = char(string(options.CubePattern));
    options.Method = lower(char(string(options.Method)));
    options.OutputPrefix = char(string(options.OutputPrefix));
    options.WriteFiles = logical(options.WriteFiles);
    options.Overwrite = logical(options.Overwrite);

    if strcmp(source.kind, 'cube')
        source.cube_pattern = options.CubePattern;
    end
end

function orbitals = load_orbitals(source)
    switch source.kind
        case 'cube'
            cube_files = dir(source.cube_pattern);
            if isempty(cube_files)
                error('calculate_scdm_spread:NoCubeFiles', ...
                    'No Cube files matching "%s" were found in %s.', ...
                    source.cube_pattern, pwd);
            end
            [~, order] = sort({cube_files.name});
            cube_files = cube_files(order);
            orbitals = repmat(empty_orbital(), numel(cube_files), 1);
            for i = 1:numel(cube_files)
                orbitals(i) = read_cube_orbital(cube_files(i).name);
            end

        case 'direct'
            orbitals = direct_orbitals(source.wan_r, source.grid_size, ...
                source.structure);

        otherwise
            error('calculate_scdm_spread:BadSource', ...
                'Unsupported source kind: %s.', source.kind);
    end
end

function orbitals = direct_orbitals(wan_r, grid_size, structure)
    validateattributes(grid_size, {'numeric'}, ...
        {'vector', 'numel', 3, 'integer', 'positive'});
    grid_size = double(grid_size(:).');
    n_grid = prod(grid_size);
    validateattributes(wan_r, {'numeric'}, {'2d', 'nonempty'});
    if size(wan_r, 1) ~= n_grid && size(wan_r, 1) ~= 2 * n_grid
        error('calculate_scdm_spread:BadWannierSize', ...
            ['wan_r has %d rows, but grid %d x %d x %d requires ', ...
             '%d scalar rows or %d spinor rows.'], ...
            size(wan_r, 1), grid_size(1), grid_size(2), grid_size(3), ...
            n_grid, 2 * n_grid);
    end

    structure = validate_direct_structure(structure);
    step_vectors = structure.supercell ./ grid_size(:);
    voxel_volume = abs(det(structure.supercell)) / n_grid;
    has_spinor = size(wan_r, 1) == 2 * n_grid;
    n_orb = size(wan_r, 2);
    orbitals = repmat(empty_orbital(), n_orb, 1);

    for iorb = 1:n_orb
        up = wan_r(1:n_grid, iorb);
        if has_spinor
            down = wan_r(n_grid + 1:end, iorb);
            amplitude_values = sqrt(abs(up).^2 + abs(down).^2) / ...
                sqrt(voxel_volume);
        else
            amplitude_values = up / sqrt(voxel_volume);
        end
        amplitude_values = reshape(amplitude_values, grid_size);
        orbitals(iorb).filename = sprintf('wan_r:%05d', iorb);
        orbitals(iorb).values = amplitude_values;
        orbitals(iorb).grid_size = grid_size;
        orbitals(iorb).origin = [0, 0, 0];
        orbitals(iorb).step_vectors = step_vectors;
        orbitals(iorb).cell = structure.supercell;
    end
end

function structure = validate_direct_structure(structure)
    if ~isstruct(structure) || ~isfield(structure, 'supercell')
        error('calculate_scdm_spread:MissingSupercell', ...
            'Direct mode structure must contain a supercell field.');
    end
    structure.supercell = double(structure.supercell);
    if ~isequal(size(structure.supercell), [3, 3]) || ...
            any(~isfinite(structure.supercell(:)))
        error('calculate_scdm_spread:BadSupercell', ...
            'structure.supercell must be a finite 3 x 3 matrix in Bohr.');
    end
end

function orbital = read_cube_orbital(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('calculate_scdm_spread:CannotOpenCube', ...
            'Cannot open Cube file %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    fgetl(fid);
    fgetl(fid);
    atom_line = sscanf(fgetl(fid), '%f');
    if numel(atom_line) < 4
        error('calculate_scdm_spread:BadCubeHeader', ...
            'Cube file %s has an invalid atom/origin line.', filename);
    end
    natoms = abs(round(atom_line(1)));
    origin = atom_line(2:4).';

    grid_size = zeros(1, 3);
    step_vectors = zeros(3, 3);
    for axis_index = 1:3
        line_values = sscanf(fgetl(fid), '%f');
        if numel(line_values) < 4
            error('calculate_scdm_spread:BadCubeGrid', ...
                'Cube file %s has an invalid grid line.', filename);
        end
        grid_size(axis_index) = abs(round(line_values(1)));
        step_vectors(axis_index, :) = line_values(2:4).';
    end

    for atom_index = 1:natoms
        atom_values = sscanf(fgetl(fid), '%f');
        if numel(atom_values) < 5
            error('calculate_scdm_spread:BadCubeAtoms', ...
                'Cube file %s has an invalid atom line.', filename);
        end
    end

    expected_values = prod(grid_size);
    raw_values = fscanf(fid, '%f');
    if numel(raw_values) < expected_values
        error('calculate_scdm_spread:BadCubeDataLength', ...
            'Cube file %s contains %d values; expected at least %d.', ...
            filename, numel(raw_values), expected_values);
    elseif numel(raw_values) > expected_values
        warning('calculate_scdm_spread:ExtraCubeData', ...
            'Cube file %s contains extra values; using the first data block.', ...
            filename);
        raw_values = raw_values(1:expected_values);
    end

    nz = grid_size(3);
    ny = grid_size(2);
    nx = grid_size(1);
    values_zyx = reshape(raw_values(1:expected_values), [nz, ny, nx]);
    values_xyz = permute(values_zyx, [3, 2, 1]);

    orbital = empty_orbital();
    orbital.filename = filename;
    orbital.values = values_xyz;
    orbital.grid_size = grid_size;
    orbital.origin = origin;
    orbital.step_vectors = step_vectors;
    orbital.cell = step_vectors .* grid_size(:);
end

function row = compute_one_spread(orbital, method)
    values = orbital.values;
    density = abs(values).^2;
    grid_size = orbital.grid_size;
    cell_matrix = orbital.cell;
    voxel_volume = abs(det(cell_matrix)) / prod(grid_size);

    weights = density(:) * voxel_volume;
    norm_value = sum(weights);
    if norm_value <= 0 || ~isfinite(norm_value)
        error('calculate_scdm_spread:BadOrbitalNorm', ...
            'Orbital %s has zero or nonfinite norm.', orbital.filename);
    end

    [ix, iy, iz] = ndgrid(0:grid_size(1)-1, ...
        0:grid_size(2)-1, 0:grid_size(3)-1);
    frac = [ix(:) / grid_size(1), iy(:) / grid_size(2), ...
        iz(:) / grid_size(3)];
    cart = orbital.origin + frac * cell_matrix;

    switch method
        case 'periodic'
            center_frac = periodic_density_center(frac, weights);
            delta_frac = frac - center_frac;
            delta_frac = delta_frac - round(delta_frac);
            delta_cart = delta_frac * cell_matrix;
            center_bohr = orbital.origin + center_frac * cell_matrix;
            spread_bohr2 = sum(weights .* sum(delta_cart.^2, 2)) / norm_value;

        case 'direct'
            center_bohr = sum(cart .* weights, 1) / norm_value;
            delta_cart = cart - center_bohr;
            spread_bohr2 = sum(weights .* sum(delta_cart.^2, 2)) / norm_value;
            center_frac = (center_bohr - orbital.origin) / cell_matrix;
            center_frac = center_frac - floor(center_frac);

        otherwise
            error('calculate_scdm_spread:BadMethod', ...
                'Unsupported spread method: %s.', method);
    end

    row = empty_spread_row();
    row.filename = orbital.filename;
    row.norm_value = norm_value;
    row.spread_bohr2 = spread_bohr2;
    row.center_frac = center_frac;
    row.center_bohr = center_bohr;
    row.max_abs = max(abs(values(:)));
    row.grid_points = prod(grid_size);
end

function center_frac = periodic_density_center(frac, weights)
    center_frac = zeros(1, 3);
    for axis_index = 1:3
        phase = sum(weights .* exp(1i * 2 * pi * frac(:, axis_index)));
        if abs(phase) <= 100 * eps(sum(weights))
            center_frac(axis_index) = 0;
        else
            center_frac(axis_index) = angle(phase) / (2 * pi);
            center_frac(axis_index) = mod(center_frac(axis_index), 1);
        end
    end
end

function [text_file, csv_file] = write_spread_outputs(spread_table, summary, options)
    text_file = [options.OutputPrefix, '.txt'];
    csv_file = [options.OutputPrefix, '.csv'];
    if ~options.Overwrite
        text_file = next_available_filename(text_file);
        csv_file = next_available_filename(csv_file);
    end

    writetable(spread_table, csv_file);

    fid = fopen(text_file, 'w');
    if fid == -1
        error('calculate_scdm_spread:CannotOpenOutput', ...
            'Cannot create output file %s.', text_file);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, 'SCDM Wannier spread report\n');
    fprintf(fid, 'Generated: %s\n', ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, 'Source: %s\n', summary.source);
    fprintf(fid, 'Method: %s\n', summary.method);
    fprintf(fid, 'Number of orbitals: %d\n', summary.n_orbitals);
    fprintf(fid, 'Total spread: %.12f Bohr^2 = %.12f Ang^2\n', ...
        summary.total_spread_bohr2, summary.total_spread_ang2);
    fprintf(fid, 'Mean spread: %.12f Bohr^2 = %.12f Ang^2\n', ...
        summary.mean_spread_bohr2, summary.mean_spread_ang2);
    fprintf(fid, 'Min/Max spread: %.12f / %.12f Ang^2\n\n', ...
        summary.min_spread_ang2, summary.max_spread_ang2);

    fprintf(fid, ['%6s  %-28s  %14s  %14s  %14s  ', ...
        '%11s %11s %11s\n'], ...
        'orb', 'file', 'norm', 'spread_A2', 'rms_A', ...
        'cf1', 'cf2', 'cf3');
    for i = 1:height(spread_table)
        fprintf(fid, ['%6d  %-28s  %14.8e  %14.8f  %14.8f  ', ...
            '%11.6f %11.6f %11.6f\n'], ...
            spread_table.orbital(i), char(spread_table.file(i)), ...
            spread_table.norm(i), spread_table.spread_ang2(i), ...
            spread_table.rms_radius_ang(i), ...
            spread_table.center_f1(i), spread_table.center_f2(i), ...
            spread_table.center_f3(i));
    end

    clear cleanup;
end

function filename = next_available_filename(filename)
    if ~isfile(filename)
        return;
    end
    [folder, stem, ext] = fileparts(filename);
    if isempty(folder)
        folder = pwd;
    end
    for index = 1:9999
        candidate = fullfile(folder, sprintf('%s_%03d%s', stem, index, ext));
        if ~isfile(candidate)
            filename = candidate;
            return;
        end
    end
    error('calculate_scdm_spread:NoAvailableFilename', ...
        'Could not find an unused output filename for %s.', filename);
end

function print_summary(summary)
    fprintf('SCDM spread calculated for %d orbitals using %s moments.\n', ...
        summary.n_orbitals, summary.method);
    fprintf('Total spread: %.8f Ang^2; mean spread: %.8f Ang^2.\n', ...
        summary.total_spread_ang2, summary.mean_spread_ang2);
    if strlength(summary.text_file) > 0
        fprintf('Spread report written to %s\n', summary.text_file);
        fprintf('Spread CSV written to %s\n', summary.csv_file);
    end
end

function orbital = empty_orbital()
    orbital = struct('filename', '', 'values', [], 'grid_size', [], ...
        'origin', [], 'step_vectors', [], 'cell', []);
end

function row = empty_spread_row()
    row = struct('filename', '', 'norm_value', NaN, 'spread_bohr2', NaN, ...
        'center_frac', [NaN, NaN, NaN], 'center_bohr', [NaN, NaN, NaN], ...
        'max_abs', NaN, 'grid_points', NaN);
end
