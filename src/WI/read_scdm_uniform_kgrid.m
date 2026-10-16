function [kpts, dims, linear_indices, grid_shift] = read_scdm_uniform_kgrid(filename, nk)
%READ_SCDM_UNIFORM_KGRID Read a complete uniform grid in reciprocal fractions.
%
% New SCDM inputs contain an explicit KGRID entry. For older files, this
% function can recover fractional coordinates from CELL_BOHR when k-points
% were accidentally stored in Cartesian reciprocal coordinates.
% grid_shift is the mesh origin in reciprocal fractional coordinates:
% kpts = integer_coords ./ dims + grid_shift (modulo reciprocal vectors).
% Its canonical representative satisfies -1/2 <= grid_shift .* dims < 1/2.
% Retain this fourth output when Fourier-interpolating a shifted mesh.

    fid = fopen(filename, 'r');
    if fid == -1
        error('read_scdm_uniform_kgrid:CannotOpenFile', ...
            'Cannot open SCDM input file %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    fgetl(fid);
    fgetl(fid);
    stored_kpts = zeros(nk, 3);
    for ik = 1:nk
        line = fgetl(fid);
        values = sscanf(line, '%g');
        if numel(values) ~= 3
            error('read_scdm_uniform_kgrid:BadKpointLine', ...
                'Invalid k-point line %d in %s.', ik, filename);
        end
        stored_kpts(ik, :) = values;
    end
    trailing_text = fread(fid, '*char').';
    clear cleanup;

    explicit_dims = parse_kgrid_metadata(trailing_text);
    supercell = parse_cell_metadata(trailing_text);

    candidates = {stored_kpts};
    candidate_names = {'stored fractional coordinates'};
    if ~isempty(supercell)
        reciprocal_cell = 2 * pi * inv(supercell).';
        candidates{end + 1} = stored_kpts / reciprocal_cell;
        candidate_names{end + 1} = 'Cartesian coordinates converted using CELL_BOHR';

        % Compatibility with older write_input_files versions, which used
        % k_cart * supercell/(2*pi) instead of the reciprocal-basis inverse.
        legacy_transform = inv(supercell).' * supercell;
        candidates{end + 1} = stored_kpts / legacy_transform;
        candidate_names{end + 1} = 'legacy nonorthogonal coordinates corrected using CELL_BOHR';
    end

    best_score = inf;
    best_name = '';
    kpts = [];
    dims = [];
    linear_indices = [];
    grid_shift = [];
    for candidate_index = 1:numel(candidates)
        candidate = candidates{candidate_index};
        if isempty(explicit_dims)
            candidate_dims = infer_grid_dims(candidate, nk);
            if isempty(candidate_dims)
                continue;
            end
        else
            candidate_dims = explicit_dims;
        end

        [valid, candidate_indices, score, candidate_shift] = map_uniform_grid( ...
            candidate, candidate_dims, nk);
        if valid && score < best_score
            best_score = score;
            best_name = candidate_names{candidate_index};
            kpts = candidate;
            dims = candidate_dims;
            linear_indices = candidate_indices;
            grid_shift = candidate_shift;
        end
    end

    if isempty(kpts)
        if isempty(explicit_dims)
            grid_description = 'no valid grid dimensions could be inferred';
        else
            grid_description = sprintf('declared KGRID is %d x %d x %d', explicit_dims);
        end
        error('read_scdm_uniform_kgrid:InvalidUniformKGrid', ...
            ['The %d k-points in %s do not form a complete uniform grid; %s. ', ...
             'Regenerate scdm.in with write_input_files.'], ...
            nk, filename, grid_description);
    end

    kpts(abs(kpts) < 1e-12) = 0;
    fprintf('  - Uniform k-grid: %d x %d x %d (%s)\n', dims, best_name);
    fprintf('  - Mesh origin (reciprocal fractions): %.12g %.12g %.12g\n', grid_shift);
end

function dims = parse_kgrid_metadata(text)
    token = regexp(text, ...
        '(?m)^\s*KGRID\s+(\d+)\s+(\d+)\s+(\d+)\s*$', ...
        'tokens', 'once');
    if isempty(token)
        dims = [];
    else
        dims = cellfun(@str2double, token);
    end
end

function supercell = parse_cell_metadata(text)
    token = regexp(text, ...
        'CELL_BOHR\s*\r?\n([^\r\n]+)\r?\n([^\r\n]+)\r?\n([^\r\n]+)', ...
        'tokens', 'once');
    if isempty(token)
        supercell = [];
        return;
    end

    supercell = zeros(3, 3);
    for row = 1:3
        values = sscanf(token{row}, '%g');
        if numel(values) ~= 3
            supercell = [];
            return;
        end
        supercell(row, :) = values;
    end
end

function dims = infer_grid_dims(kpts, nk)
    rounded = round(kpts * 1e8) / 1e8;
    dims = zeros(1, 3);
    for axis = 1:3
        dims(axis) = numel(unique(rounded(:, axis)));
    end
    if prod(dims) ~= nk
        dims = [];
    end
end

function [valid, linear_indices, residual, grid_shift] = map_uniform_grid(kpts, dims, nk)
    valid = false;
    linear_indices = [];
    residual = inf;
    grid_shift = [];
    if numel(dims) ~= 3 || prod(dims) ~= nk || any(dims < 1)
        return;
    end

    scaled = kpts .* dims;
    % Use the same half-step representative for positive/negative MP points.
    % This keeps origin and indexing stable under input reordering.
    offset = mod(scaled(1, :) + 0.5, 1) - 0.5;
    offset(abs(offset) < 1e-10) = 0;
    offset(abs(abs(offset) - 0.5) < 1e-10) = -0.5;
    integer_coords = round(scaled - offset);
    residual = max(abs(scaled - offset - integer_coords), [], 'all');
    if residual > 1e-6
        return;
    end

    indices = mod(integer_coords, dims) + 1;
    linear_indices = sub2ind(dims, indices(:, 1), indices(:, 2), indices(:, 3));
    valid = numel(unique(linear_indices)) == nk;
    grid_shift = offset ./ dims;
end
