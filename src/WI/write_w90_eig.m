function filename = write_w90_eig(eig_matrix, seedname, conv_factor)
%WRITE_W90_EIG Write KSSOLV eigenvalues to a Wannier90 .eig file.
%
% filename = write_w90_eig(eig_matrix, seedname)
% filename = write_w90_eig(eig_matrix, seedname, conv_factor)
%
% eig_matrix must have size [nbnd, nkpts]. Energies written to the .eig
% file must be in eV. Use conv_factor = 27.211386 for Hartree input or
% conv_factor = 1 when eig_matrix is already in eV.

if nargin < 2
    error('write_w90_eig:NotEnoughInputs', ...
        'eig_matrix and seedname are required.');
end

validateattributes(eig_matrix, {'numeric'}, ...
    {'2d', 'nonempty', 'real', 'finite'}, ...
    'write_w90_eig', 'eig_matrix');
if ~(ischar(seedname) || (isstring(seedname) && isscalar(seedname)))
    error('write_w90_eig:BadSeedname', ...
        'seedname must be a character vector or string scalar.');
end

if nargin < 3 || isempty(conv_factor)
    % Unit inference is inherently ambiguous; explicit conv_factor is safer.
    if max(abs(eig_matrix(:))) <= 5
        conv_factor = 27.211386;
        warning('write_w90_eig:InferredHartree', ...
            ['Input energies were inferred to be Hartree and will be ', ...
             'converted to eV. Pass conv_factor explicitly to override.']);
    else
        conv_factor = 1;
        warning('write_w90_eig:InferredEV', ...
            ['Input energies were inferred to be eV. Pass conv_factor ', ...
             'explicitly to override.']);
    end
else
    validateattributes(conv_factor, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        'write_w90_eig', 'conv_factor');
end

seedname = char(seedname);
if length(seedname) < 4 || ~strcmpi(seedname(end-3:end), '.eig')
    filename = [seedname, '.eig'];
else
    filename = seedname;
end

output_dir = fileparts(filename);
if ~isempty(output_dir) && ~exist(output_dir, 'dir')
    [ok, message] = mkdir(output_dir);
    if ~ok
        error('write_w90_eig:CannotCreateDirectory', ...
            'Cannot create output directory: %s', message);
    end
end

[fid, message] = fopen(filename, 'w');
if fid < 0
    error('write_w90_eig:CannotOpenFile', ...
        'Cannot open %s for writing: %s', filename, message);
end
cleanup = onCleanup(@() fclose(fid));

[nbnd, nkpts] = size(eig_matrix);

% Wannier90 expects all bands of k-point 1, then all bands of k-point 2.
for k = 1:nkpts
    for b = 1:nbnd
        energy_ev = eig_matrix(b, k) * conv_factor;
        fprintf(fid, '%5d %5d %18.12f\n', b, k, energy_ev);
    end
end

clear cleanup
end
