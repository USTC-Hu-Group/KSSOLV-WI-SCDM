function [filename, info] = write_w90_mmn(seedname, varargin)
%WRITE_W90_MMN Export a Wannier90 .mmn file from KSSOLV UNK files.
%
% [filename, info] = write_w90_mmn(seedname)
% [filename, info] = write_w90_mmn(seedname, 'SCDMInput', 'scdm.in', ...)
%
% The overlap matrix is written in Wannier90 order:
%   ik ikp G1 G2 G3
%   M(1,1), M(2,1), ..., M(nbnd,1), M(1,2), ...
%
% For an nnkp entry k+b = kp+G, this routine evaluates
%   M_mn = <u_mk | exp(-i G.r) | u_n,kp>
% on the same real-space grid and normalization convention used by scdm_m.

if nargin < 1 || isempty(seedname)
    seedname = 'wannier90';
end
seedname = char(seedname);

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'SCDMInput', 'scdm.in', @(x) ischar(x) || isstring(x));
addParameter(parser, 'NNKPFile', [seedname, '.nnkp'], @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputFile', [seedname, '.mmn'], @(x) ischar(x) || isstring(x));
addParameter(parser, 'UnkDirectory', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'NumBands', [], @(x) isempty(x) || isnumeric(x));
addParameter(parser, 'CacheSize', 12, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(parser, 'NormalizeUNK', true, @(x) islogical(x) || isnumeric(x));
parse(parser, varargin{:});

scdm_input = char(parser.Results.SCDMInput);
nnkp_file = char(parser.Results.NNKPFile);
filename = char(parser.Results.OutputFile);
unk_dir = char(parser.Results.UnkDirectory);
num_bands_override = parser.Results.NumBands;
cache_size = max(0, round(parser.Results.CacheSize));
normalize_unk = logical(parser.Results.NormalizeUNK);

if ~isfile(scdm_input)
    error('write_w90_mmn:MissingSCDMInput', ...
        'Cannot find SCDM input file: %s', scdm_input);
end
if ~isfile(nnkp_file)
    error('write_w90_mmn:MissingNNKP', ...
        'Cannot find Wannier90 nnkp file: %s', nnkp_file);
end

[input_dir, ~, ~] = fileparts(scdm_input);
if isempty(input_dir)
    input_dir = pwd;
end
if isempty(unk_dir)
    unk_dir = fullfile(input_dir, 'tmp');
end

[nband, nk] = read_scdm_counts(scdm_input);
if ~isempty(num_bands_override)
    validateattributes(num_bands_override, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
        mfilename, 'NumBands');
    nband = double(num_bands_override);
end

nnkp = read_w90_nnkp(nnkp_file, nk);
ensure_parent_dir(filename);

cache_keys = nan(cache_size, 1);
cache_values = cell(cache_size, 1);
cache_dims = cell(cache_size, 1);
cache_count = 0;
next_cache_slot = 1;

fid = fopen(filename, 'w');
if fid < 0
    error('write_w90_mmn:CannotOpenOutput', ...
        'Cannot open MMN output file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Created by write_w90_mmn on %s\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, '%10d %10d %10d\n', nband, nk, nnkp.nntot);

timer = tic;
for ik = 1:nk
    [psi_k, dims_k] = get_unk_cached(ik);
    for ib = 1:nnkp.nntot
        ikp = nnkp.nnlist(ik, ib);
        gvec = squeeze(nnkp.nncell(ik, ib, :)).';
        [psi_kp, dims_kp] = get_unk_cached(ikp);
        if ~isequal(dims_k, dims_kp)
            error('write_w90_mmn:GridMismatch', ...
                'UNK grid mismatch between k=%d and k=%d.', ik, ikp);
        end

        mmn = compute_mmn_block(psi_k, psi_kp, gvec, dims_k);
        fprintf(fid, '%5d %5d %5d %5d %5d\n', ...
            ik, ikp, gvec(1), gvec(2), gvec(3));
        write_mmn_block(fid, mmn);
    end
    if mod(ik, 25) == 0 || ik == nk
        fprintf('MMN exported for %d / %d k-points\n', ik, nk);
    end
end

elapsed = toc(timer);
clear cleanup;

info = struct();
info.filename = filename;
info.sc_dm_input = scdm_input;
info.nnkp_file = nnkp_file;
info.unk_dir = unk_dir;
info.num_bands = nband;
info.nkpts = nk;
info.nntot = nnkp.nntot;
info.elapsed_seconds = elapsed;
fprintf('MMN file written: %s (%.3f seconds)\n', filename, elapsed);

    function [psi, dims] = get_unk_cached(ik_query)
        hit = find(cache_keys(1:cache_count) == ik_query, 1);
        if ~isempty(hit)
            psi = cache_values{hit};
            dims = cache_dims{hit};
            return;
        end

        [psi, dims] = read_w90_unk_spinor(ik_query, nband, input_dir, ...
            unk_dir, normalize_unk);

        if cache_size > 0
            if cache_count < cache_size
                cache_count = cache_count + 1;
                cache_keys(cache_count) = ik_query;
                cache_values{cache_count} = psi;
                cache_dims{cache_count} = dims;
            else
                cache_keys(next_cache_slot) = ik_query;
                cache_values{next_cache_slot} = psi;
                cache_dims{next_cache_slot} = dims;
                next_cache_slot = next_cache_slot + 1;
                if next_cache_slot > cache_size
                    next_cache_slot = 1;
                end
            end
        end
    end
end

function [nband, nk] = read_scdm_counts(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('write_w90_mmn:CannotOpenSCDMInput', ...
        'Cannot open SCDM input file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));
values = sscanf(fgetl(fid), '%f');
if numel(values) < 2
    error('write_w90_mmn:InvalidSCDMHeader', ...
        'The first line of %s must contain at least num_bands and nkpts.', filename);
end
nband = round(values(1));
nk = round(values(2));
end

function nnkp = read_w90_nnkp(filename, nk)
lines = regexp(fileread(filename), '\r\n|\n|\r', 'split');
begin_idx = find(strcmpi(strtrim(lines), 'begin nnkpts'), 1);
if isempty(begin_idx)
    error('write_w90_mmn:MissingNNKPTSBlock', ...
        'Cannot find begin nnkpts block in %s.', filename);
end

line_idx = begin_idx + 1;
[line_idx, line] = next_data_line(lines, line_idx);
nntot_values = sscanf(line, '%d');
if numel(nntot_values) ~= 1 || nntot_values(1) <= 0
    error('write_w90_mmn:InvalidNNTOT', ...
        'Invalid nntot value in %s.', filename);
end
nntot = double(nntot_values(1));

nnlist = zeros(nk, nntot);
nncell = zeros(nk, nntot, 3);
for ik = 1:nk
    for ib = 1:nntot
        [line_idx, line] = next_data_line(lines, line_idx + 1);
        values = sscanf(line, '%d');
        if numel(values) < 5
            error('write_w90_mmn:InvalidNNKPEntry', ...
                'Invalid nnkpts entry in %s near line %d.', filename, line_idx);
        end
        if values(1) ~= ik
            warning('write_w90_mmn:UnexpectedNNKPIndex', ...
                'nnkp entry line %d starts with %d, expected %d. Using Wannier90 sequential order.', ...
                line_idx, values(1), ik);
        end
        ikp = values(2);
        if ikp < 1 || ikp > nk
            error('write_w90_mmn:InvalidNeighborIndex', ...
                'Invalid neighbor k-point index %d in %s.', ikp, filename);
        end
        nnlist(ik, ib) = ikp;
        nncell(ik, ib, :) = values(3:5);
    end
end

nnkp = struct();
nnkp.nntot = nntot;
nnkp.nnlist = nnlist;
nnkp.nncell = nncell;
end

function [line_idx, line] = next_data_line(lines, line_idx)
while line_idx <= numel(lines)
    line = strtrim(lines{line_idx});
    if ~isempty(line) && ~startsWith(line, '!') && ~startsWith(line, '#')
        if startsWith(lower(line), 'end ')
            error('write_w90_mmn:UnexpectedEndBlock', ...
                'Unexpected end block while reading nnkpts.');
        end
        return;
    end
    line_idx = line_idx + 1;
end
error('write_w90_mmn:UnexpectedEOF', ...
    'Unexpected end of file while reading nnkpts.');
end

function [unk, dims] = read_w90_unk_spinor(kpoint_index, nband_expected, input_dir, ...
        unk_dir, normalize_unk)
filename = find_unk_file(kpoint_index, input_dir, unk_dir);
fid = fopen(filename, 'r');
if fid < 0
    error('write_w90_mmn:CannotOpenUNK', ...
        'Cannot open UNK file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

nx = fscanf(fid, '%d', 1);
ny = fscanf(fid, '%d', 1);
nz = fscanf(fid, '%d', 1);
fscanf(fid, '%d', 1);
nband_file = fscanf(fid, '%d', 1);

if isempty(nx) || isempty(ny) || isempty(nz) || isempty(nband_file)
    error('write_w90_mmn:InvalidUNKHeader', ...
        'Invalid UNK header in %s.', filename);
end
if nband_expected > nband_file
    error('write_w90_mmn:TooFewUNKBands', ...
        'UNK file %s contains %d bands, but %d are required.', ...
        filename, nband_file, nband_expected);
end

n_grid = nx * ny * nz;
raw_data = fscanf(fid, '%g');
if numel(raw_data) == 4 * n_grid * nband_file
    psi_up = complex(raw_data(1:4:end), raw_data(2:4:end));
    psi_dw = complex(raw_data(3:4:end), raw_data(4:4:end));
elseif numel(raw_data) == 2 * n_grid * nband_file
    psi_up = complex(raw_data(1:2:end), raw_data(2:2:end));
    psi_dw = zeros(size(psi_up));
else
    error('write_w90_mmn:InvalidUNKDataSize', ...
        'File %s has %d data values, expected %d for SOC or %d for scalar UNK.', ...
        filename, numel(raw_data), 4 * n_grid * nband_file, ...
        2 * n_grid * nband_file);
end

psi_up = reshape(psi_up, [n_grid, nband_file]);
psi_dw = reshape(psi_dw, [n_grid, nband_file]);
psi_up = psi_up(:, 1:nband_expected);
psi_dw = psi_dw(:, 1:nband_expected);

if normalize_unk
    for band = 1:nband_expected
        norm_factor = norm([psi_up(:, band); psi_dw(:, band)]);
        if norm_factor > 1e-12
            psi_up(:, band) = psi_up(:, band) / norm_factor;
            psi_dw(:, band) = psi_dw(:, band) / norm_factor;
        end
    end
end

unk = [psi_up; psi_dw];
dims = [nx, ny, nz];
end

function filename = find_unk_file(kpoint_index, input_dir, unk_dir)
base = sprintf('UNK%05d', kpoint_index);
candidates = {
    fullfile(unk_dir, [base '.NC'])
    fullfile(unk_dir, [base '.1'])
    fullfile(input_dir, 'tmp', [base '.NC'])
    fullfile(input_dir, 'tmp', [base '.1'])
    fullfile(input_dir, [base '.NC'])
    fullfile(input_dir, [base '.1'])
    fullfile(input_dir, 'UNK_files', [base '.NC'])
    fullfile(input_dir, 'UNK_files', [base '.1'])
    fullfile(pwd, 'tmp', [base '.NC'])
    fullfile(pwd, 'tmp', [base '.1'])
    fullfile(pwd, [base '.NC'])
    fullfile(pwd, [base '.1'])
};

for i = 1:numel(candidates)
    if isfile(candidates{i})
        filename = candidates{i};
        return;
    end
end
error('write_w90_mmn:MissingUNK', ...
    'Cannot find UNK file for k-point %d. Expected %s.NC or %s.1.', ...
    kpoint_index, base, base);
end

function mmn = compute_mmn_block(psi_k, psi_kp, gvec, dims)
if size(psi_k, 1) ~= size(psi_kp, 1) || size(psi_k, 2) ~= size(psi_kp, 2)
    error('write_w90_mmn:WavefunctionSizeMismatch', ...
        'Wavefunction matrices for neighboring k-points have different sizes.');
end

if all(gvec == 0)
    mmn = psi_k' * psi_kp;
    return;
end

nx = dims(1);
ny = dims(2);
nz = dims(3);
[x_grid, y_grid, z_grid] = ndgrid(0:nx-1, 0:ny-1, 0:nz-1);
r_frac = [x_grid(:) / nx, y_grid(:) / ny, z_grid(:) / nz];
phase = exp(-1i * 2 * pi * (r_frac * double(gvec(:))));
phase = [phase; phase];
mmn = psi_k' * (phase .* psi_kp);
end

function write_mmn_block(fid, mmn)
nband = size(mmn, 1);
for n = 1:nband
    for m = 1:nband
        fprintf(fid, '%18.12f %18.12f\n', real(mmn(m, n)), imag(mmn(m, n)));
    end
end
end

function ensure_parent_dir(filename)
parent = fileparts(filename);
if ~isempty(parent) && ~isfolder(parent)
    mkdir(parent);
end
end
