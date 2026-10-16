function cube_files = write_wannier_cubes(wan_r, grid_size, structure, seedname, output_dir)
%WRITE_WANNIER_CUBES Write real-space Wannier orbitals as Gaussian Cube files.
%
% Scalar Wannier functions are phase-aligned and written as real-valued
% amplitudes. For spinor Wannier functions, one scalar Cube file cannot
% represent both complex components, so the total spinor amplitude
% sqrt(|up|^2 + |down|^2) is written.

    if nargin < 5 || isempty(output_dir)
        output_dir = pwd;
    end
    if nargin < 4 || isempty(seedname)
        seedname = 'scdm';
    end

    validateattributes(grid_size, {'numeric'}, ...
        {'vector', 'numel', 3, 'integer', 'positive'});
    grid_size = double(grid_size(:).');
    nx = grid_size(1);
    ny = grid_size(2);
    nz = grid_size(3);
    n_grid = prod(grid_size);

    validateattributes(wan_r, {'numeric'}, {'2d', 'nonempty'});
    if size(wan_r, 1) ~= n_grid && size(wan_r, 1) ~= 2 * n_grid
        error('write_wannier_cubes:BadWannierSize', ...
            ['wan_r has %d rows, but grid %d x %d x %d requires ', ...
             '%d scalar rows or %d spinor rows.'], ...
            size(wan_r, 1), nx, ny, nz, n_grid, 2 * n_grid);
    end

    structure = validate_structure(structure);
    seedname = sanitize_seedname(seedname);
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    voxel_volume = abs(det(structure.supercell)) / n_grid;
    amplitude_scale = 1 / sqrt(voxel_volume);
    n_wannier = size(wan_r, 2);
    cube_files = cell(n_wannier, 1);
    has_spinor_rows = size(wan_r, 1) == 2 * n_grid;

    for iw = 1:n_wannier
        up = wan_r(1:n_grid, iw);
        if has_spinor_rows
            down = wan_r(n_grid + 1:end, iw);
        else
            down = zeros(size(up));
        end

        down_fraction = norm(down) / max(norm([up; down]), eps);
        if down_fraction > 1e-10
            values = sqrt(abs(up).^2 + abs(down).^2);
            description = 'total spinor amplitude sqrt(|up|^2+|down|^2)';
        else
            [values, imaginary_fraction] = phase_aligned_real_part(up);
            description = 'phase-aligned real Wannier amplitude';
            if imaginary_fraction > 1e-5
                warning('write_wannier_cubes:ComplexScalarOrbital', ...
                    ['Wannier function %d retains an imaginary norm fraction ', ...
                     'of %.3e after phase alignment; only its real part is written.'], ...
                    iw, imaginary_fraction);
            end
        end

        values = values * amplitude_scale;
        filename = sprintf('%s_%05d.cube', seedname, iw);
        cube_files{iw} = fullfile(output_dir, filename);
        write_one_cube(cube_files{iw}, values, grid_size, structure, ...
            seedname, iw, description);
    end

    fprintf('Gaussian Cube output written: %d Wannier orbitals, seedname %s\n', ...
        n_wannier, seedname);
end

function structure = validate_structure(structure)
    required = {'supercell', 'atomic_numbers', 'effective_charges', 'xyzlist'};
    for index = 1:numel(required)
        if ~isfield(structure, required{index})
            error('write_wannier_cubes:MissingStructureField', ...
                'Structure metadata is missing field %s.', required{index});
        end
    end

    structure.supercell = double(structure.supercell);
    structure.atomic_numbers = double(structure.atomic_numbers(:));
    structure.effective_charges = double(structure.effective_charges(:));
    structure.xyzlist = double(structure.xyzlist);
    natoms = numel(structure.atomic_numbers);

    if ~isequal(size(structure.supercell), [3, 3]) || ...
            ~isequal(size(structure.xyzlist), [natoms, 3]) || ...
            numel(structure.effective_charges) ~= natoms
        error('write_wannier_cubes:BadStructureSize', ...
            'Invalid supercell, atom-number, charge, or coordinate dimensions.');
    end
    if any(~isfinite(structure.supercell(:))) || ...
            any(~isfinite(structure.xyzlist(:))) || ...
            any(~isfinite(structure.effective_charges)) || ...
            any(~isfinite(structure.atomic_numbers))
        error('write_wannier_cubes:NonfiniteStructure', ...
            'Structure metadata contains NaN or Inf values.');
    end
end

function seedname = sanitize_seedname(seedname)
    seedname = char(string(seedname));
    [~, seedname, ~] = fileparts(seedname);
    seedname = regexprep(seedname, '[^A-Za-z0-9_.-]', '_');
    if isempty(seedname)
        seedname = 'scdm';
    end
end

function [values, imaginary_fraction] = phase_aligned_real_part(values)
    phase_measure = sum(values.^2);
    if abs(phase_measure) > eps
        values = values * exp(-0.5i * angle(phase_measure));
    end
    imaginary_fraction = norm(imag(values)) / max(norm(values), eps);
    values = real(values);
end

function write_one_cube(filename, values, grid_size, structure, seedname, ...
        wannier_index, description)

    fid = fopen(filename, 'w');
    if fid == -1
        error('write_wannier_cubes:CannotOpenFile', ...
            'Cannot create Cube file %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    natoms = numel(structure.atomic_numbers);
    step_vectors = structure.supercell ./ grid_size(:);

    fprintf(fid, '%s Wannier function %05d\n', seedname, wannier_index);
    fprintf(fid, 'Generated by MATLAB SCDM; %s; units: Bohr\n', description);
    fprintf(fid, '%5d %16.10f %16.10f %16.10f\n', natoms, 0, 0, 0);
    for axis = 1:3
        fprintf(fid, '%5d %16.10f %16.10f %16.10f\n', ...
            grid_size(axis), step_vectors(axis, :));
    end

    for atom_index = 1:natoms
        fprintf(fid, '%5d %16.10f %16.10f %16.10f %16.10f\n', ...
            structure.atomic_numbers(atom_index), ...
            structure.effective_charges(atom_index), ...
            structure.xyzlist(atom_index, :));
    end

    volume = reshape(values, grid_size);
    z_fastest_values = permute(volume, [3, 2, 1]);
    z_fastest_values = z_fastest_values(:);
    for first = 1:6:numel(z_fastest_values)
        last = min(first + 5, numel(z_fastest_values));
        fprintf(fid, ' %13.5E', z_fastest_values(first:last));
        fprintf(fid, '\n');
    end

    clear cleanup;
end
