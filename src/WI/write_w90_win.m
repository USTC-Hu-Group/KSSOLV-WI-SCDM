function filename = write_w90_win(sys, num_wann, seedname, varargin)
%WRITE_W90_WIN Write a Wannier90 .win file from a KSSOLV Crystal object.
%
% filename = write_w90_win(sys, num_wann, seedname)
% filename = write_w90_win(..., 'LengthUnit', 'ang')
% filename = write_w90_win(..., 'AtomCoordinates', 'frac')
% filename = write_w90_win(..., 'KpointCoordinates', 'frac')
% filename = write_w90_win(..., 'LatticeVectors', 'columns')
%
% KSSOLV conventions used by the default "auto" modes:
%   * supercell lattice vectors are rows and lengths are in Bohr.
%   * xyzlist contains Cartesian coordinates.
%   * kpts of a finalized Crystal are Cartesian reciprocal coordinates.
%
% For a plain struct, "auto" assumes fractional k points and identifies
% atom coordinates in [0,1] as fractional. Override the relevant option
% when this assumption is not appropriate.

if nargin < 2
    error('write_w90_win:NotEnoughInputs', ...
        'sys and num_wann are required.');
end

if nargin < 3 || isempty(seedname)
    seedname = read_member(sys, {'name'}, false);
    if isempty(seedname)
        seedname = 'wannier90';
    end
    seedname = regexprep(char(seedname), '[^A-Za-z0-9_.-]', '_');
end

parser = inputParser;
parser.FunctionName = 'write_w90_win';
addParameter(parser, 'LengthUnit', 'bohr', @is_text_scalar);
addParameter(parser, 'AtomCoordinates', 'auto', @is_text_scalar);
addParameter(parser, 'KpointCoordinates', 'auto', @is_text_scalar);
addParameter(parser, 'LatticeVectors', 'rows', @is_text_scalar);
parse(parser, varargin{:});

length_unit = lower(char(parser.Results.LengthUnit));
atom_mode = lower(char(parser.Results.AtomCoordinates));
kpoint_mode = lower(char(parser.Results.KpointCoordinates));
lattice_mode = lower(char(parser.Results.LatticeVectors));

must_be_one_of(length_unit, {'bohr', 'ang'}, 'LengthUnit');
must_be_one_of(atom_mode, {'auto', 'cart', 'frac'}, 'AtomCoordinates');
must_be_one_of(kpoint_mode, {'auto', 'cart', 'frac'}, 'KpointCoordinates');
must_be_one_of(lattice_mode, {'rows', 'columns'}, 'LatticeVectors');

validateattributes(num_wann, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    'write_w90_win', 'num_wann');

lattice = double(read_member(sys, {'supercell'}, true));
if ~isequal(size(lattice), [3, 3])
    error('write_w90_win:BadSupercell', ...
        'sys.supercell must be a 3-by-3 matrix.');
end
if strcmp(lattice_mode, 'columns')
    lattice = lattice.';
end
if abs(det(lattice)) < 100 * eps(max(abs(lattice(:))))
    error('write_w90_win:SingularSupercell', ...
        'The lattice-vector matrix is singular.');
end

nbnd = double(read_member(sys, {'nbnd'}, true));
validateattributes(nbnd, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    'write_w90_win', 'sys.nbnd');
if num_wann > nbnd
    error('write_w90_win:TooManyWannierFunctions', ...
        'num_wann (%d) cannot exceed num_bands (%d).', num_wann, nbnd);
end

% Noncollinear Bloch states are singly occupied two-component spinors.
% Leave scalar output unchanged; SOC calculations must explicitly tell W90.
nspin = read_member(sys, {'nspin'}, false);
noncolin = read_member(sys, {'noncolin'}, false);
is_spinor = isequal(nspin, 4) || isequal(noncolin, true) || isequal(noncolin, 1);

xyzlist = double(read_member(sys, {'xyzlist'}, true));
if size(xyzlist, 2) ~= 3
    error('write_w90_win:BadCoordinates', ...
        'sys.xyzlist must have exactly three columns.');
end
natom = size(xyzlist, 1);
atom_labels = get_atom_labels(sys, natom);

% atoms_frac avoids any ambiguity in the length unit of atomic positions.
if strcmp(atom_mode, 'auto')
    if isa(sys, 'Crystal') || isa(sys, 'Molecule')
        atom_mode = 'cart';
    elseif coordinates_look_fractional(xyzlist)
        atom_mode = 'frac';
    else
        atom_mode = 'cart';
    end
end
if strcmp(atom_mode, 'cart')
    atom_frac = xyzlist / lattice;
else
    atom_frac = xyzlist;
end
atom_frac(abs(atom_frac) < 1e-13) = 0;

kpoints = double(read_member(sys, {'kpts', 'k-points'}, true));
if size(kpoints, 2) < 3
    error('write_w90_win:BadKpoints', ...
        'The k-point list must have at least three columns.');
end
kpoints = kpoints(:, 1:3); % Explicitly discard optional k-point weights.

if strcmp(kpoint_mode, 'auto')
    if isa(sys, 'Crystal')
        kpoint_mode = 'cart';
    else
        kpoint_mode = 'frac';
    end
end
if strcmp(kpoint_mode, 'cart')
    % finalize(Crystal) uses k_cart = k_frac * 2*pi*inv(supercell)'.
    kpoints_frac = kpoints * lattice.' / (2 * pi);
else
    kpoints_frac = kpoints;
end
kpoints_frac(abs(kpoints_frac) < 1e-13) = 0;

mp_grid = determine_mp_grid(sys, kpoints_frac);

lattice_out = lattice;
if strcmp(length_unit, 'ang')
    lattice_out = lattice_out * 0.529177210903;
end

seedname = char(seedname);
if length(seedname) < 4 || ~strcmpi(seedname(end-3:end), '.win')
    filename = [seedname, '.win'];
else
    filename = seedname;
end

output_dir = fileparts(filename);
if ~isempty(output_dir) && ~exist(output_dir, 'dir')
    [ok, message] = mkdir(output_dir);
    if ~ok
        error('write_w90_win:CannotCreateDirectory', ...
            'Cannot create output directory: %s', message);
    end
end

[fid, message] = fopen(filename, 'w');
if fid < 0
    error('write_w90_win:CannotOpenFile', ...
        'Cannot open %s for writing: %s', filename, message);
end
cleanup = onCleanup(@() fclose(fid));

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
fprintf(fid, '! Generated by write_w90_win on %s\n', timestamp);
fprintf(fid, '! Source: KSSOLV crystal data\n\n');
fprintf(fid, 'num_wann  = %d\n', num_wann);
fprintf(fid, 'num_bands = %d\n\n', nbnd);
if is_spinor
    fprintf(fid, 'spinors = true\n\n');
end

fprintf(fid, 'begin unit_cell_cart\n');
fprintf(fid, '%s\n', length_unit);
write_three_columns(fid, lattice_out);
fprintf(fid, 'end unit_cell_cart\n\n');

fprintf(fid, 'begin atoms_frac\n');
for iatom = 1:natom
    fprintf(fid, '%-6s %20.12f %20.12f %20.12f\n', ...
        atom_labels{iatom}, atom_frac(iatom, 1), ...
        atom_frac(iatom, 2), atom_frac(iatom, 3));
end
fprintf(fid, 'end atoms_frac\n\n');

fprintf(fid, 'mp_grid = %d %d %d\n\n', mp_grid(1), mp_grid(2), mp_grid(3));
fprintf(fid, 'begin kpoints\n');
write_three_columns(fid, kpoints_frac);
fprintf(fid, 'end kpoints\n');

clear cleanup
end

function value = read_member(object, names, required)
value = [];
for iname = 1:numel(names)
    name = names{iname};
    try
        if isstruct(object) && isfield(object, name)
            value = object.(name);
        elseif isobject(object) && isprop(object, name)
            value = object.(name);
        else
            value = [];
        end
    catch
        value = [];
    end
    if ~isempty(value)
        return;
    end

    try
        value = get(object, name);
    catch
        value = [];
    end
    if ~isempty(value)
        return;
    end
end

if required
    error('write_w90_win:MissingField', ...
        'Cannot find required sys property: %s.', strjoin(names, ' or '));
end
end

function labels = get_atom_labels(sys, natom)
atomlist = read_member(sys, {'atomlist'}, false);
if ~isempty(atomlist)
    labels = labels_from_atom_array(atomlist, natom);
    return;
end

atoms = read_member(sys, {'atoms'}, false);
alist = read_member(sys, {'alist'}, false);
if isempty(atoms) || isempty(alist)
    error('write_w90_win:MissingAtoms', ...
        'Cannot determine atom labels from atomlist or atoms/alist.');
end
alist = double(alist(:));
if numel(alist) ~= natom || any(alist < 1) || any(alist ~= round(alist))
    error('write_w90_win:BadAtomIndex', ...
        'sys.alist must contain one valid species index per atom.');
end

labels = cell(natom, 1);
for iatom = 1:natom
    labels{iatom} = atom_symbol(index_array(atoms, alist(iatom)));
end
end

function labels = labels_from_atom_array(atomlist, natom)
if ischar(atomlist)
    labels = cellstr(atomlist);
elseif isstring(atomlist)
    labels = cellstr(atomlist(:));
else
    if numel(atomlist) ~= natom
        error('write_w90_win:BadAtomList', ...
            'atomlist must contain one entry per row of xyzlist.');
    end
    labels = cell(natom, 1);
    for iatom = 1:natom
        labels{iatom} = atom_symbol(index_array(atomlist, iatom));
    end
end

if numel(labels) ~= natom
    error('write_w90_win:BadAtomList', ...
        'atomlist must contain one entry per row of xyzlist.');
end
labels = labels(:);
end

function item = index_array(array, index)
if iscell(array)
    item = array{index};
else
    item = array(index);
end
end

function symbol = atom_symbol(atom)
if ischar(atom)
    symbol = strtrim(atom);
elseif isstring(atom) && isscalar(atom)
    symbol = char(atom);
else
    try
        symbol = char(atom.symbol);
    catch
        error('write_w90_win:MissingAtomSymbol', ...
            'Each atom/species entry must provide a symbol.');
    end
end
if isempty(symbol)
    error('write_w90_win:EmptyAtomSymbol', 'An atom symbol is empty.');
end
end

function mp_grid = determine_mp_grid(sys, kpoints_frac)
mp_grid = read_member(sys, {'nkxyz'}, false);
if ~isempty(mp_grid)
    mp_grid = double(mp_grid(:).');
end

if ~valid_mp_grid(mp_grid, size(kpoints_frac, 1))
    tol = 1e-10;
    rounded = round(kpoints_frac / tol) * tol;
    inferred = [numel(unique(rounded(:, 1))), ...
        numel(unique(rounded(:, 2))), numel(unique(rounded(:, 3)))];
    if valid_mp_grid(inferred, size(kpoints_frac, 1))
        mp_grid = inferred;
    else
        n123 = [read_member(sys, {'n1'}, false), ...
            read_member(sys, {'n2'}, false), ...
            read_member(sys, {'n3'}, false)];
        if valid_mp_grid(n123, size(kpoints_frac, 1))
            mp_grid = double(n123);
        else
            error('write_w90_win:CannotDetermineMPGrid', ...
                ['Cannot determine a valid Monkhorst-Pack grid. Set ', ...
                 'sys.nkxyz or provide a complete regular k-point mesh.']);
        end
    end
end
end

function tf = valid_mp_grid(grid, nkpoint)
tf = isnumeric(grid) && numel(grid) == 3 && ...
    all(isfinite(grid)) && all(grid > 0) && ...
    all(grid == round(grid)) && prod(grid) == nkpoint;
end

function tf = coordinates_look_fractional(coordinates)
tol = 1e-10;
tf = all(coordinates(:) >= -tol) && all(coordinates(:) <= 1 + tol);
end

function write_three_columns(fid, values)
for irow = 1:size(values, 1)
    fprintf(fid, '%20.12f %20.12f %20.12f\n', ...
        values(irow, 1), values(irow, 2), values(irow, 3));
end
end

function must_be_one_of(value, choices, option_name)
if ~any(strcmp(value, choices))
    error('write_w90_win:BadOption', ...
        '%s must be one of: %s.', option_name, strjoin(choices, ', '));
end
end

function tf = is_text_scalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
