function [wan_r, Amn, cube_files] = scdm_m(fin, fout_amn)
%SCDM_M SCDM projections from cell-periodic scalar/spinor UNK states.
%
% Inputs
%   fin      Input file containing band count, k-points, eigenvalue file,
%            named filtering mode, and optional chemical potential/smearing.
%   fout_amn Optional output filename for the Wannier90 .amn file.
%
% Outputs
%   wan_r    Real-space Wannier functions, stacked as [spin up; spin down].
%   Amn      Cell array of projection matrices, one [nband, nwan] matrix per k.
%   cube_files Cell array containing generated Gaussian Cube filenames.
%
% The SOC UNK format is assumed to be text with interleaved values:
%   Re_up Im_up Re_down Im_down
% UNK contains u_nk, with psi_nk(r) = exp(i*k.r) u_nk(r).
% The anchor is the input point minimizing sum(abs(k_fractional)); ties
% use input order. No separate Gamma calculation or manual anchor is used.

scdm_timer = tic;

[input_dir, ~, ~] = fileparts(fin);
if isempty(input_dir)
    input_dir = pwd;
end

if nargin < 2
    fout_amn = 'scdm.amn';
end
cube_files = {};

[nband, nk, nwan_input, feig, focc, mu, sigma, kpts, structure] = ...
    read_scdm_input(fin);
if ~isfile(feig)
    feig = fullfile(input_dir, feig);
end

E = load(feig);
eigKpts = reshape(E(:,3), nband, nk);
occ = compute_occupations(eigKpts, focc, mu, sigma);
report_filter_diagnostics(occ, focc, mu, sigma, nband, nwan_input);

[~, anchor_idx] = min(sum(abs(kpts), 2));
if any(kpts(anchor_idx, :) ~= 0)
    warning('scdm_m:NonGammaAnchor', ...
        ['Gamma point not strictly found. Using k-point %d with ', ...
         'minimum fractional-coordinate L1 norm.'], anchor_idx);
end
fprintf('SCDM anchor: index %d, reciprocal fractions [%.12g %.12g %.12g].\n', ...
    anchor_idx, kpts(anchor_idx, :));

[Psi_anchor, nx, ny, nz] = readUNK_spinor(anchor_idx, nband, input_dir);
n_grid = nx * ny * nz;

disp('Computing SCDM columns (QRCP)...');
weighted_anchor = Psi_anchor .* reshape(occ(:, anchor_idx), 1, []);
[~, R, perm] = qr(weighted_anchor', 'vector');
report_weighted_rank(R, nwan_input, focc, anchor_idx);
cols = select_columns_from_qr(R, perm, nwan_input);
nwan = numel(cols);

[Xs, Ys, Zs] = ndgrid(0:nx-1, 0:ny-1, 0:nz-1);
rposFrac = [Xs(:) / nx, Ys(:) / ny, Zs(:) / nz];

disp('Performing projections and constructing Wannier functions...');
Amn = cell(nk, 1);
wan_r_total = zeros(2 * n_grid, nwan);

for k = 1:nk
    spatial_cols = mod(cols - 1, n_grid) + 1;
    centerFrac = rposFrac(spatial_cols, :);
    centerFrac = centerFrac - round(centerFrac);
    phase_factor_proj = exp(1i * 2 * pi * centerFrac * kpts(k, :)');

    Psi = readUNK_spinor(k, nband, input_dir);
    Psi_proj = diag(phase_factor_proj) * Psi(cols, :) * diag(occ(:, k));

    [Umn, ~, Vmn] = svd(Psi_proj', 0);
    Amn{k} = Umn * Vmn';

    Psi_mixed = Psi * Amn{k};
    phase_grid = exp(1i * 2 * pi * (rposFrac * kpts(k, :)'));
    full_phase_grid = [phase_grid; phase_grid];
    wan_r_total = wan_r_total + Psi_mixed .* full_phase_grid;

    if mod(k, 50) == 0 || k == nk
        fprintf('Processed %d / %d k-points\n', k, nk);
    end
end

wan_r = wan_r_total / nk;

write_amn_soc(Amn, nk, fout_amn);
[amn_dir, seedname, ~] = fileparts(fout_amn);
if isempty(amn_dir)
    amn_dir = pwd;
end
if isempty(seedname)
    seedname = 'scdm';
end

if isempty(structure)
    warning('scdm_m:MissingCubeMetadata', ...
        ['No Cube structure metadata was found in %s. Re-run ', ...
         'write_input_files to enable Gaussian Cube output.'], fin);
    disp('Calculation done. AMN file written; Cube output skipped.');
else
    cube_files = write_wannier_cubes(wan_r, [nx, ny, nz], structure, ...
        seedname, amn_dir);
    disp('Calculation done. AMN and Gaussian Cube outputs written.');
end

scdm_runtime_seconds = toc(scdm_timer);
runtime_file = fullfile(input_dir, 'scdm_runtime.txt');
write_scdm_runtime(runtime_file, scdm_runtime_seconds, fin, fout_amn);
fprintf('scdm_m runtime: %.3f seconds (%s).\n', ...
    scdm_runtime_seconds, format_elapsed_time(scdm_runtime_seconds));

end

function [nband, nk, nwan, feig, focc, mu, sigma, kpts, structure] = ...
        read_scdm_input(fin)
fid = fopen(fin, 'r');
if fid == -1
    error('Cannot open input file %s', fin);
end
cleanup = onCleanup(@() fclose(fid));

header = sscanf(fgetl(fid), '%g');
if numel(header) ~= 3
    error('First line of %s must contain nband, nk, and nwan.', fin);
end
nband = header(1);
nk = header(2);
nwan = header(3);

filter_line = strsplit(strtrim(fgetl(fid)));
if numel(filter_line) < 2
    error('Second line of %s must contain eigenvalue file and filtering mode.', fin);
end
feig = filter_line{1};
focc = normalize_filter_mode(filter_line{2});

if strcmp(focc, 'isolated')
    mu = [];
    sigma = [];
else
    if numel(filter_line) ~= 4
        error('%s mode requires mu and sigma on the second line of %s.', ...
            focc, fin);
    end
    mu = str2double(filter_line{3});
    sigma = str2double(filter_line{4});
    if ~isfinite(mu) || ~isfinite(sigma) || sigma <= 0
        error('Invalid mu or sigma on the second line of %s.', fin);
    end
end

kpts = zeros(nk, 3);
for k = 1:nk
    kpts(k, :) = fscanf(fid, '%g', 3);
end
fgetl(fid);
structure = read_cube_metadata(fid, fin);
[kpts, ~, ~] = read_scdm_uniform_kgrid(fin, nk);
end

function structure = read_cube_metadata(fid, fin)
    structure = [];
    line = next_nonempty_line(fid);
    if ~ischar(line)
        return;
    end
    if ~strcmp(strtrim(line), 'BEGIN_CUBE_METADATA')
        warning('scdm_m:UnknownTrailingInput', ...
            'Ignoring unrecognized trailing content in %s.', fin);
        return;
    end
    structure = struct();

    line = next_nonempty_line(fid);
    if ischar(line) && startsWith(strtrim(line), 'KGRID ')
        kgrid_values = sscanf(line, 'KGRID %g %g %g');
        if numel(kgrid_values) ~= 3
            error('Invalid KGRID metadata in %s.', fin);
        end
        structure.kgrid = kgrid_values(:).';
        line = next_nonempty_line(fid);
    end
    if ~ischar(line) || ~strcmp(strtrim(line), 'CELL_BOHR')
        error('Cube metadata in %s is missing CELL_BOHR.', fin);
    end

    supercell = zeros(3, 3);
    for row = 1:3
        values = sscanf(fgetl(fid), '%g');
        if numel(values) ~= 3
            error('Invalid CELL_BOHR row %d in %s.', row, fin);
        end
        supercell(row, :) = values;
    end

    atom_line = next_nonempty_line(fid);
    if ~ischar(atom_line)
        error('Cube metadata in %s is missing ATOMS_BOHR.', fin);
    end
    atom_header = strsplit(strtrim(atom_line));
    if numel(atom_header) ~= 2 || ~strcmp(atom_header{1}, 'ATOMS_BOHR')
        error('Cube metadata in %s is missing ATOMS_BOHR.', fin);
    end
    natoms = str2double(atom_header{2});
    if ~isfinite(natoms) || natoms < 0 || natoms ~= round(natoms)
        error('Invalid atom count in Cube metadata of %s.', fin);
    end

    atomic_numbers = zeros(natoms, 1);
    effective_charges = zeros(natoms, 1);
    xyzlist = zeros(natoms, 3);
    for atom_index = 1:natoms
        values = sscanf(fgetl(fid), '%g');
        if numel(values) ~= 5
            error('Invalid atom row %d in Cube metadata of %s.', atom_index, fin);
        end
        atomic_numbers(atom_index) = values(1);
        effective_charges(atom_index) = values(2);
        xyzlist(atom_index, :) = values(3:5);
    end

    line = next_nonempty_line(fid);
    if ~ischar(line) || ~strcmp(strtrim(line), 'END_CUBE_METADATA')
        error('Cube metadata in %s is missing END_CUBE_METADATA.', fin);
    end

    structure.supercell = supercell;
    structure.atomic_numbers = atomic_numbers;
    structure.effective_charges = effective_charges;
    structure.xyzlist = xyzlist;
end

function line = next_nonempty_line(fid)
    line = fgetl(fid);
    while ischar(line) && isempty(strtrim(line))
        line = fgetl(fid);
    end
end

function occ = compute_occupations(eigKpts, focc, mu, sigma)
switch focc
    case 'gaussian'
        occ = exp(-((eigKpts - mu).^2) / (sigma^2));
    case 'erfc'
        occ = 0.5 * erfc((eigKpts - mu) / sigma);
    case 'isolated'
        occ = ones(size(eigKpts));
    otherwise
        error('Unsupported SCDM filtering mode: %s.', focc);
end
end

function report_filter_diagnostics(occ, focc, mu, sigma, nband, nwan_input)
if strcmp(focc, 'isolated')
    fprintf('SCDM filter: isolated; all band weights are 1.\n');
else
    fprintf(['SCDM filter: %s, mu = %.8f eV, sigma = %.8f eV; ', ...
        'weight range = [%.3e, %.3e].\n'], ...
        focc, mu, sigma, min(occ(:)), max(occ(:)));
end

if ~strcmp(focc, 'isolated') && nwan_input >= nband
    warning('scdm_m:FullBandSpaceFilter', ...
        ['num_wann=%d equals or exceeds num_bands=%d. The %s filter is ', ...
         'applied and can change the SCDM gauge/Wannier functions, but the ', ...
         'retained band subspace is complete. Interpolated band eigenvalues ', ...
         'are therefore expected to be nearly insensitive to mu and sigma.'], ...
        round(nwan_input), nband, focc);
end
end

function report_weighted_rank(R, nwan_input, focc, anchor_idx)
diagR = abs(diag(R));
if isempty(diagR)
    effective_rank = 0;
else
    tolerance = max(size(R)) * eps(max(diagR));
    effective_rank = sum(diagR > max(tolerance, 1e-12 * max(diagR)));
end
fprintf('Weighted SCDM matrix effective rank at k=%d: %d.\n', ...
    anchor_idx, effective_rank);

if nwan_input >= 1 && effective_rank < round(nwan_input)
    warning('scdm_m:RankDeficientFilter', ...
        ['The %s-weighted SCDM matrix has effective rank %d, below the ', ...
         'requested num_wann=%d. Some selected directions are controlled ', ...
         'by numerically tiny weights; adjust mu/sigma or reduce num_wann.'], ...
        focc, effective_rank, round(nwan_input));
end
end

function write_scdm_runtime(filename, elapsed_seconds, fin, fout_amn)
fid = fopen(filename, 'w');
if fid == -1
    warning('scdm_m:CannotWriteRuntime', ...
        'Cannot write scdm_m runtime record to %s.', filename);
    return;
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'timestamp %s\n', char(datetime('now', ...
    'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'elapsed_seconds %.12f\n', elapsed_seconds);
fprintf(fid, 'elapsed_hms %s\n', format_elapsed_time(elapsed_seconds));
fprintf(fid, 'input_file %s\n', fin);
fprintf(fid, 'amn_file %s\n', fout_amn);
clear cleanup;
end

function text = format_elapsed_time(elapsed_seconds)
hours = floor(elapsed_seconds / 3600);
minutes = floor(mod(elapsed_seconds, 3600) / 60);
seconds = mod(elapsed_seconds, 60);
text = sprintf('%02d:%02d:%06.3f', hours, minutes, seconds);
end

function mode = normalize_filter_mode(value)
mode = lower(strtrim(char(value)));
valid_modes = {'isolated', 'erfc', 'gaussian'};
if ~any(strcmp(mode, valid_modes))
    error(['SCDM filtering mode must be isolated, erfc, or gaussian; ', ...
        'numeric focc values are no longer supported.']);
end
end

function cols = select_columns_from_qr(R, perm, nwan_input)
if nwan_input >= 1
    nwan = round(nwan_input);
else
    diagR = abs(diag(R));
    if isempty(diagR) || diagR(1) == 0
        error('Cannot determine nwan from QR diagonal because R is singular.');
    end

    cutoff_idx = find(diagR / diagR(1) < nwan_input, 1);
    if isempty(cutoff_idx)
        nwan = numel(diagR);
    else
        nwan = max(cutoff_idx - 1, 1);
    end
end

if nwan > numel(perm)
    error('Requested %d Wannier functions, but only %d candidate columns are available.', ...
        nwan, numel(perm));
end

cols = perm(1:nwan);
end

function [Unk, nx, ny, nz] = readUNK_spinor(k, nband_expected, input_dir)
fname = find_unk_file(k, input_dir);
fid = fopen(fname, 'r');
if fid == -1
    error('Cannot open file %s', fname);
end
cleanup = onCleanup(@() fclose(fid));

nx = fscanf(fid, '%d', 1);
ny = fscanf(fid, '%d', 1);
nz = fscanf(fid, '%d', 1);
fscanf(fid, '%d', 1);
nband = fscanf(fid, '%d', 1);

if nband ~= nband_expected
    warning('Expected %d bands, but %s header has %d bands.', ...
        nband_expected, fname, nband);
end

n_grid = nx * ny * nz;
raw_data = fscanf(fid, '%g');

if numel(raw_data) == 4 * n_grid * nband
    psi_up = complex(raw_data(1:4:end), raw_data(2:4:end));
    psi_dw = complex(raw_data(3:4:end), raw_data(4:4:end));
elseif numel(raw_data) == 2 * n_grid * nband
    psi_up = complex(raw_data(1:2:end), raw_data(2:2:end));
    psi_dw = zeros(size(psi_up));
else
    error('File %s has %d data values, expected %d for SOC or %d for scalar UNK.', ...
        fname, numel(raw_data), 4 * n_grid * nband, 2 * n_grid * nband);
end

psi_up = reshape(psi_up, [n_grid, nband]);
psi_dw = reshape(psi_dw, [n_grid, nband]);

for band = 1:nband
    norm_factor = sqrt(norm(psi_up(:, band))^2 + norm(psi_dw(:, band))^2);
    if norm_factor > 1e-10
        psi_up(:, band) = psi_up(:, band) / norm_factor;
        psi_dw(:, band) = psi_dw(:, band) / norm_factor;
    end
end

Unk = [psi_up; psi_dw];
end

function fname = find_unk_file(k, input_dir)
base = sprintf('UNK%05d', k);
candidates = {
    fullfile(input_dir, 'tmp', [base '.NC'])
    fullfile(input_dir, 'tmp', [base '.1'])
    fullfile('tmp', [base '.NC'])
    fullfile('tmp', [base '.1'])
    fullfile(input_dir, [base '.NC'])
    fullfile(input_dir, 'UNK_files', [base '.NC'])
    fullfile(input_dir, [base '.1'])
    fullfile(input_dir, 'UNK_files', [base '.1'])
    [base '.NC']
    fullfile('UNK_files', [base '.NC'])
    [base '.1']
    fullfile('UNK_files', [base '.1'])
};

for idx = 1:numel(candidates)
    if isfile(candidates{idx})
        fname = candidates{idx};
        return;
    end
end

error('Cannot open UNK file for k-point %d. Tried %s/.NC and %s/.1 in tmp, %s, and UNK_files.', ...
    k, base, base, input_dir);
end

function write_amn_soc(Amn, nk, fname)
fid = fopen(fname, 'w');
if fid == -1
    error('Cannot open AMN output file %s', fname);
end
cleanup = onCleanup(@() fclose(fid));

[nb, nw] = size(Amn{1});
fprintf(fid, 'Created by MATLAB SCDM-SOC\n');
fprintf(fid, '%10d %10d %10d\n', nb, nk, nw);

for k = 1:nk
    Ak = Amn{k};
    for w = 1:nw
        for b = 1:nb
            fprintf(fid, '%5d %5d %5d %16.12f %16.12f\n', ...
                b, w, k, real(Ak(b, w)), imag(Ak(b, w)));
        end
    end
end
end
