function [proj_vals, pao_info] = compute_pao_projectability(crystal, nbnd, nk, count_only, bloch_basis)
%COMPUTE_PAO_PROJECTABILITY Compute Bloch-state projectability onto PAOs.
%
% proj_vals(n,k) = <psi_nk|P_A(k)|psi_nk>, where P_A(k) projects onto
% the orthonormalized pseudo-atomic-orbital (PAO) subspace at k.
% pao_info.total_count is the raw PAO projection dimension N_PAO. For SOC
% calculations it includes both spinor components.
%
% PAO radial functions are read from PP_PSWFC/PP_CHI.* in the ppdata
% sublibrary selected by kssolvpptype. If available, crystal.ppvar.funct is
% used as the functional actually loaded from the pseudopotential headers;
% crystal.funct is only a fallback calculation label.
% Bloch wavefunctions are read from ./tmp/UNKxxxxx.1 or UNKxxxxx.NC.
% An optional fifth input accepts a KSSOLV Basis object and reads its
% reciprocal coefficients directly, avoiding the large UNK text round trip.

if nargin < 4
    count_only = false;
end
if nargin < 5
    bloch_basis = [];
end

unk_dir = fullfile(pwd, 'tmp');

if ~count_only && isempty(bloch_basis) && ~isfolder(unk_dir)
    error('compute_pao_projectability:MissingUNKDirectory', ...
        '未找到 UNK 文件目录：%s。请先运行 write_UNK_files。', unk_dir);
end

validateattributes(nbnd, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(nk, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});

[symbols, ~, xyzlist] = crystal_atom_data(crystal);
funct = char(read_crystal_member(crystal, 'funct', true));
[pao_funct, pao_funct_source] = infer_pseudopotential_functional( ...
    crystal, funct);
nspin = read_crystal_member(crystal, 'nspin', false);
is_soc = ~isempty(nspin) && double(nspin) == 4;

unique_symbols = unique(symbols, 'stable');
[pp_root, pp_library, pp_type, pp_extension, selected_pp_files] = resolve_ppdata_root( ...
    unique_symbols, pao_funct, is_soc);
fprintf('PAO pseudopotential root selected from MATLAB path: %s\n', pp_root);
fprintf('PAO pseudopotential library locked by kssolvpptype("%s","%s"): %s\n', ...
    pp_type, pp_extension, pp_library);
fprintf('PAO pseudopotential functional requested by crystal.funct: %s\n', funct);
if ~strcmpi(strtrim(pao_funct), strtrim(funct))
    fprintf('PAO pseudopotential functional inferred from %s: %s\n', ...
        pao_funct_source, pao_funct);
end

pao_by_symbol = cell(numel(unique_symbols), 1);
for itype = 1:numel(unique_symbols)
    pp_file = selected_pp_files{itype};
    pao_by_symbol{itype} = read_upf_pao(pp_file);
    fprintf('PAO pseudopotential used: %s <- %s\n', ...
        unique_symbols{itype}, pp_file);
end

% Map every atom to the PAO data selected for its chemical species.
pao_by_atom = cell(numel(symbols), 1);
for iatom = 1:numel(symbols)
    type_index = find(strcmp(unique_symbols, symbols{iatom}), 1);
    pao_by_atom{iatom} = pao_by_symbol{type_index};
end

n_spatial_pao = 0;
for iatom = 1:numel(pao_by_atom)
    channels = pao_by_atom{iatom}.channels;
    for ichannel = 1:numel(channels)
        n_spatial_pao = n_spatial_pao + 2 * channels(ichannel).l + 1;
    end
end
if n_spatial_pao == 0
    error('compute_pao_projectability:NoPAO', ...
        '所选赝势中没有可用的 PP_CHI.* 伪原子轨道。');
end

pao_info.spatial_count = n_spatial_pao;
pao_info.spin_factor = 1 + double(is_soc);
pao_info.total_count = n_spatial_pao * pao_info.spin_factor;
pao_info.rank_by_k = zeros(1, nk);
pao_info.pp_root = pp_root;
pao_info.pp_library = pp_library;
pao_info.pp_type = pp_type;
pao_info.pp_extension = pp_extension;
pao_info.pp_files = selected_pp_files;
pao_info.crystal_funct = funct;
pao_info.pp_functional = pao_funct;
pao_info.pp_functional_source = pao_funct_source;
pao_info.count_only = logical(count_only);
pao_info.ranks_computed = ~logical(count_only);
pao_info.species = unique_symbols;
pao_info.atom_symbols = symbols;
pao_info.channels_by_species = pao_by_symbol;
pao_info.metric = 'identity in the retained plane-wave basis; no PAW/US augmentation';
pao_info.spinor_construction = 'spatial PAOs tensor spin; no j,mj recoupling';
pao_info.singular_values_by_k = cell(1, nk);
pao_info.rank_tolerance_by_k = nan(1, nk);
pao_info.raw_projectability = [];

proj_vals = zeros(nbnd, nk);
if count_only
    proj_vals = [];
    return;
end

[pw_indices, g_vectors] = plane_wave_basis(crystal);
pao_info.pw_indices = pw_indices;
pao_info.raw_projectability = nan(nbnd, nk);
kpoints_cart = double(read_crystal_member(crystal, 'kpts', true));
kpoints_cart = kpoints_cart(:, 1:3);
if size(kpoints_cart, 1) ~= nk
    error('compute_pao_projectability:BadKpoints', ...
        'crystal.kpts 的行数与 nk 不一致。');
end

for ik = 1:nk
    q_vectors = g_vectors + kpoints_cart(ik, :);
    pao_coeff = build_pao_coefficients( ...
        q_vectors, xyzlist, pao_by_atom, n_spatial_pao);

    if isempty(bloch_basis)
        [psi_coeff, unk_is_soc] = read_unk_reciprocal( ...
            unk_dir, ik, nbnd, pw_indices);
    else
        [psi_coeff, unk_is_soc] = read_basis_reciprocal( ...
            bloch_basis, ik, nbnd, pw_indices, crystal);
    end
    if unk_is_soc
        spatial_coeff = pao_coeff;
        pao_coeff = complex(zeros( ...
            2 * size(spatial_coeff, 1), 2 * size(spatial_coeff, 2)));
        pao_coeff(1:size(spatial_coeff, 1), ...
            1:size(spatial_coeff, 2)) = spatial_coeff;
        pao_coeff(size(spatial_coeff, 1) + 1:end, ...
            size(spatial_coeff, 2) + 1:end) = spatial_coeff;
    elseif is_soc
        error('compute_pao_projectability:SpinMismatch', ...
            'crystal.nspin=4，但 UNK%05d 不是旋量波函数。', ik);
    end

    % SVD gives an orthonormal basis Q for the PAO Hilbert subspace.
    % The projector is P_A = Q*Q', avoiding explicit inverse overlap.
    [Q, singular_values, ~] = svd(pao_coeff, 'econ');
    singular_values = diag(singular_values);
    if isempty(singular_values) || singular_values(1) <= 0
        error('compute_pao_projectability:SingularPAOSpace', ...
            'k 点 %d 的 PAO 子空间为空或奇异。', ik);
    end
    rank_tolerance = max(size(pao_coeff)) * eps(singular_values(1));
    keep = singular_values > max(rank_tolerance, 1e-10 * singular_values(1));
    pao_info.singular_values_by_k{ik} = singular_values;
    pao_info.rank_tolerance_by_k(ik) = max(rank_tolerance, 1e-10 * singular_values(1));
    Q = Q(:, keep);
    if isempty(Q)
        error('compute_pao_projectability:SingularPAOSpace', ...
            'k 点 %d 的 PAO 子空间数值秩为零。', ik);
    end
    pao_info.rank_by_k(ik) = size(Q, 2);

    psi_norm = sqrt(sum(abs(psi_coeff).^2, 1));
    if any(psi_norm < 1e-12)
        error('compute_pao_projectability:ZeroWavefunction', ...
            'k 点 %d 存在零范数 Bloch 波函数。', ik);
    end
    psi_coeff = psi_coeff ./ psi_norm;

    p = real(sum(abs(Q' * psi_coeff).^2, 1));
    pao_info.raw_projectability(:, ik) = p(:);
    proj_vals(:, ik) = min(max(p(:), 0), 1);

    if mod(ik, 100) == 0 || ik == 1 || ik == nk
        fprintf(['已计算 PAO 可投影性：k = %d/%d，PAO 子空间秩 = %d，', ...
            'p 范围 = [%.6f, %.6f]\n'], ...
            ik, nk, size(Q, 2), min(p), max(p));
    end
end
end

function [psi_coeff, is_soc] = read_basis_reciprocal( ...
        bloch_basis, ik, nbnd, indices, crystal)
% Read KSSOLV reciprocal coefficients directly, equivalently to the
% write_UNK_files -> FFT round trip, without creating large text files.
if isprop(bloch_basis, 'wavefuncell') && ~isempty(bloch_basis.wavefuncell)
    wave_data = bloch_basis.wavefuncell{ik, 1};
    psi_raw = wave_data.psi;
    source_indices = double(wave_data.idxnz(:));
elseif ik == 1 && isprop(bloch_basis, 'psi')
    psi_raw = bloch_basis.psi;
    source_indices = double(bloch_basis.idxnz(:));
else
    error('compute_pao_projectability:BadBlochBasis', ...
        'KSSOLV Bloch basis does not contain k point %d.', ik);
end
if size(psi_raw, 2) < nbnd
    error('compute_pao_projectability:TooFewBasisBands', ...
        'KSSOLV basis at k=%d contains %d bands; %d are required.', ...
        ik, size(psi_raw, 2), nbnd);
end

n_grid = double(crystal.n1) * double(crystal.n2) * double(crystal.n3);
npw = numel(source_indices);
is_soc = size(psi_raw, 1) == 2 * npw;
if ~(size(psi_raw, 1) == npw || is_soc)
    error('compute_pao_projectability:BadBasisShape', ...
        'Unexpected reciprocal coefficient shape at k=%d.', ik);
end
if any(source_indices < 1) || any(source_indices > n_grid) || ...
        any(indices < 1) || any(indices > n_grid)
    error('compute_pao_projectability:BadBasisIndices', ...
        'Plane-wave indices exceed the FFT grid at k=%d.', ik);
end

if is_soc
    up_grid = complex(zeros(n_grid, nbnd));
    down_grid = complex(zeros(n_grid, nbnd));
    up_grid(source_indices, :) = psi_raw(1:npw, 1:nbnd);
    down_grid(source_indices, :) = psi_raw(npw + 1:end, 1:nbnd);
    psi_coeff = [up_grid(indices, :); down_grid(indices, :)];
else
    reciprocal_grid = complex(zeros(n_grid, nbnd));
    reciprocal_grid(source_indices, :) = psi_raw(:, 1:nbnd);
    psi_coeff = reciprocal_grid(indices, :);
end
end

function coeff = build_pao_coefficients(q_vectors, xyzlist, pao_by_atom, n_pao)
nq = size(q_vectors, 1);
max_l = 0;
for iatom = 1:numel(pao_by_atom)
    max_l = max(max_l, max([pao_by_atom{iatom}.channels.l]));
end
y_lm = ylmr(q_vectors(:, 1), q_vectors(:, 2), q_vectors(:, 3), max_l);
q_norm = sqrt(sum(q_vectors.^2, 2));

coeff = complex(zeros(nq, n_pao));
offset = 0;
for iatom = 1:numel(pao_by_atom)
    atom_phase = exp(-1i * (q_vectors * xyzlist(iatom, :).'));
    channels = pao_by_atom{iatom}.channels;
    for ichannel = 1:numel(channels)
        channel = channels(ichannel);
        l = channel.l;
        radial_q = radial_bessel_transform( ...
            channel.r, channel.rab, channel.chi, l, q_norm);
        lm_indices = l^2 + (1:(2 * l + 1));
        columns = offset + (1:(2 * l + 1));
        coeff(:, columns) = ((-1i)^l .* radial_q .* atom_phase) ...
            .* y_lm(:, lm_indices);
        offset = offset + 2 * l + 1;
    end
end
end

function values = radial_bessel_transform(r, rab, chi, l, q)
% PP_CHI stores u_l(r)=r*R_l(r), hence the transform contains r*chi.
chunk_size = 512;
values = zeros(numel(q), 1);
radial_factor = r(:) .* chi(:);
for first = 1:chunk_size:numel(q)
    last = min(first + chunk_size - 1, numel(q));
    q_chunk = q(first:last).';
    jl = bessel(l, r(:) * q_chunk);
    integrand = radial_factor .* jl;
    values(first:last) = simpson(numel(r), integrand, rab(:)).';
end
end

function [psi_coeff, is_soc] = read_unk_reciprocal(unk_dir, ik, nbnd, indices)
nc_file = fullfile(unk_dir, sprintf('UNK%05d.NC', ik));
scalar_file = fullfile(unk_dir, sprintf('UNK%05d.1', ik));
if isfile(nc_file)
    filename = nc_file;
    is_soc = true;
elseif isfile(scalar_file)
    filename = scalar_file;
    is_soc = false;
else
    error('compute_pao_projectability:MissingUNK', ...
        '找不到 k 点 %d 对应的 UNK 文件。', ik);
end

fid = fopen(filename, 'r');
if fid < 0
    error('compute_pao_projectability:CannotOpenUNK', ...
        '无法打开 %s。', filename);
end
cleanup = onCleanup(@() fclose(fid));

nx = fscanf(fid, '%d', 1);
ny = fscanf(fid, '%d', 1);
nz = fscanf(fid, '%d', 1);
fscanf(fid, '%d', 1);
file_nbnd = fscanf(fid, '%d', 1);
if file_nbnd < nbnd
    error('compute_pao_projectability:TooFewBands', ...
        '%s 仅包含 %d 条能带，需要 %d 条。', filename, file_nbnd, nbnd);
end

n_grid = nx * ny * nz;
if any(indices < 1) || any(indices > n_grid)
    error('compute_pao_projectability:BadPWIndices', ...
        '平面波索引超出 UNK FFT 网格范围。');
end

raw = fscanf(fid, '%g');
if is_soc
    expected = 4 * n_grid * file_nbnd;
    if numel(raw) ~= expected
        error('compute_pao_projectability:BadUNKData', ...
            '%s 的旋量数据长度错误。', filename);
    end
    complex_data = complex(raw(1:2:end), raw(2:2:end));
    complex_data = reshape(complex_data, 2, n_grid, file_nbnd);
    psi_up = reshape(complex_data(1, :, 1:nbnd), n_grid, nbnd);
    psi_down = reshape(complex_data(2, :, 1:nbnd), n_grid, nbnd);
    psi_coeff = [fft_selected(psi_up, nx, ny, nz, indices); ...
        fft_selected(psi_down, nx, ny, nz, indices)];
else
    expected = 2 * n_grid * file_nbnd;
    if numel(raw) ~= expected
        error('compute_pao_projectability:BadUNKData', ...
            '%s 的标量数据长度错误。', filename);
    end
    complex_data = complex(raw(1:2:end), raw(2:2:end));
    psi_real = reshape(complex_data, n_grid, file_nbnd);
    psi_coeff = fft_selected(psi_real(:, 1:nbnd), nx, ny, nz, indices);
end

clear cleanup
end

function selected = fft_selected(psi_real, nx, ny, nz, indices)
selected = complex(zeros(numel(indices), size(psi_real, 2)));
for ibnd = 1:size(psi_real, 2)
    reciprocal_grid = fftn(reshape(psi_real(:, ibnd), nx, ny, nz));
    selected(:, ibnd) = reciprocal_grid(indices);
end
end

function [indices, g_vectors] = plane_wave_basis(crystal)
use_qe_indices = false;
try
    use_qe_indices = logical(crystal.modifyidxnz) && ~isempty(crystal.qeidxnz);
catch
end

grid = Ggrid(crystal);
if use_qe_indices
    indices = double(crystal.qeidxnz(:));
else
    indices = double(grid.idxnz(:));
end
g_vectors = double(grid.kkxyz(indices, :));
end

function pao = read_upf_pao(filename)
text = fileread(filename);
r = read_upf_numeric_block(text, 'PP_R');
rab = read_upf_numeric_block(text, 'PP_RAB');

tokens = regexp(text, ...
    '<PP_CHI\.\d+([^>]*)>([\s\S]*?)</PP_CHI\.\d+>', 'tokens');
if isempty(tokens)
    error('compute_pao_projectability:NoPAOInUPF', ...
        '赝势 %s 不包含非空 PP_CHI.* 数据。', filename);
end

channels = repmat(struct('l', 0, 'label', '', 'r', [], ...
    'rab', [], 'chi', []), numel(tokens), 1);
for ichannel = 1:numel(tokens)
    attributes = tokens{ichannel}{1};
    l_token = regexp(attributes, '\sl\s*=\s*"([^"]+)"', 'tokens', 'once');
    label_token = regexp(attributes, ...
        '\slabel\s*=\s*"([^"]+)"', 'tokens', 'once');
    if isempty(l_token)
        error('compute_pao_projectability:MissingAngularMomentum', ...
            '赝势 %s 的 PP_CHI.* 缺少角动量 l。', filename);
    end

    chi = sscanf(tokens{ichannel}{2}, '%f');
    if numel(chi) ~= numel(r)
        error('compute_pao_projectability:BadPAOSize', ...
            '赝势 %s 的 PP_CHI.* 与 PP_R 网格长度不一致。', filename);
    end

    channels(ichannel).l = str2double(l_token{1});
    if isempty(label_token)
        channels(ichannel).label = sprintf('l=%d', channels(ichannel).l);
    else
        channels(ichannel).label = strtrim(label_token{1});
    end
    channels(ichannel).r = r;
    channels(ichannel).rab = rab;
    channels(ichannel).chi = chi;
end
pao.channels = channels;
end

function values = read_upf_numeric_block(text, tag)
pattern = ['<' tag '(?:\s[^>]*)?>([\s\S]*?)</' tag '>'];
token = regexp(text, pattern, 'tokens', 'once');
if isempty(token)
    error('compute_pao_projectability:MissingUPFBlock', ...
        'UPF 文件缺少 %s 数据块。', tag);
end
values = sscanf(token{1}, '%f');
end

function [root, library, pp_type, pp_extension, selected_files] = ...
        resolve_ppdata_root(symbols, funct, is_soc)
%RESOLVE_PPDATA_ROOT Lock a ppdata library using the current kssolvpptype.

    candidates = find_ppdata_candidates_on_path();
    if isempty(candidates)
        error('compute_pao_projectability:MissingPpDirectory', ...
            ['No directory named ppdata could be found from the MATLAB path. ', ...
             'Add the KSSOLV tree or its ppdata directory with addpath first.']);
    end

    if exist('kssolvpptype', 'file') ~= 2
        error('compute_pao_projectability:MissingKssolvPpType', ...
            ['kssolvpptype is not on the MATLAB path. Run KSSOLV_startup ', ...
             'before computing projectability.']);
    end
    [pp_type, pp_extension] = kssolvpptype();
    pp_type = strtrim(char(pp_type));
    pp_extension = lower(strtrim(char(pp_extension)));
    if ~strcmpi(pp_extension, 'upf')
        error('compute_pao_projectability:UnsupportedPpExtension', ...
            ['PAO projectability requires UPF PP_CHI.* data, but the ', ...
             'current kssolvpptype extension is "%s".'], pp_extension);
    end

    library_name = pp_type;
    if strcmpi(pp_type, 'default')
        if is_soc
            library_name = 'default_FR';
        else
            library_name = 'default';
        end
    end

    failure_messages = cell(numel(candidates), 1);
    for candidate_index = 1:numel(candidates)
        candidate = candidates{candidate_index};
        files = cell(numel(symbols), 1);
        try
            candidate_library = find_pptype_subfolder( ...
                candidate, library_name, pp_type, pp_extension, symbols);
            for symbol_index = 1:numel(symbols)
                files{symbol_index} = find_matching_pao_upf( ...
                    candidate_library, symbols{symbol_index}, funct, is_soc);
            end
            root = candidate;
            library = candidate_library;
            selected_files = files;
            return;
        catch candidate_error
            failure_messages{candidate_index} = sprintf('%s: %s', ...
                candidate, candidate_error.message);
        end
    end

    error('compute_pao_projectability:NoUsablePpDirectory', ...
        ['Found %d ppdata directories from the MATLAB path, but none ', ...
         'contained the kssolvpptype-selected library "%s" with matching ', ...
         'PAO UPF files for all species.\n%s'], ...
        numel(candidates), library_name, strjoin(failure_messages, newline));
end

function library = find_pptype_subfolder(root, library_name, pp_type, ...
        pp_extension, symbols)
%FIND_PPTYPE_SUBFOLDER Map kssolvpptype to one immediate ppdata sublibrary.

    library_name = strtrim(char(library_name));
    if isempty(library_name) || contains(library_name, filesep) || ...
            (~ispc && contains(library_name, '\')) || contains(library_name, '/')
        error('compute_pao_projectability:BadPpTypeFolder', ...
            'kssolvpptype must select a nonempty ppdata subfolder name.');
    end

    entries = dir(root);
    entries = entries([entries.isdir]);
    names = {entries.name};
    keep = ~ismember(names, {'.', '..'});
    entries = entries(keep);
    names = names(keep);

    match = find(strcmpi(names, library_name));
    if isempty(match) && ~strcmpi(pp_type, 'default')
        inferred_matches = false(size(entries));
        for entry_index = 1:numel(entries)
            folder = fullfile(root, entries(entry_index).name);
            inferred_matches(entry_index) = folder_provides_pptype( ...
                folder, symbols, pp_type, pp_extension);
        end
        match = find(inferred_matches);
    end

    if isempty(match)
        if isempty(names)
            available = '(none)';
        else
            available = strjoin(sort(names), ', ');
        end
        error('compute_pao_projectability:MissingPpTypeFolder', ...
            ['ppdata root %s does not contain the subfolder requested by ', ...
             'kssolvpptype="%s". Available subfolders: %s'], ...
            root, library_name, available);
    end
    if numel(match) > 1
        error('compute_pao_projectability:AmbiguousPpTypeFolder', ...
            ['ppdata root %s contains multiple subfolders matching ', ...
             'kssolvpptype="%s" case-insensitively.'], root, library_name);
    end

    library = absolute_folder(fullfile(root, entries(match).name));
end

function tf = folder_provides_pptype(folder, symbols, pp_type, pp_extension)
%FOLDER_PROVIDES_PPTYPE Require matching files for every chemical species.

    files = dir(folder);
    files = files(~[files.isdir]);
    filenames = {files.name};
    tf = true;
    for symbol_index = 1:numel(symbols)
        expression = ['^' regexptranslate('escape', symbols{symbol_index}) ...
            '[._-]' regexptranslate('escape', pp_type) ...
            '[\w.-]*\.' regexptranslate('escape', pp_extension) '$'];
        if ~any(~cellfun(@isempty, regexpi(filenames, expression, 'once')))
            tf = false;
            return;
        end
    end
end

function candidates = find_ppdata_candidates_on_path()
%FIND_PPDATA_CANDIDATES_ON_PATH Search path entries and their ancestors.

    path_entries = [{pwd}, strsplit(path, pathsep)];
    candidates = {};
    candidate_keys = {};

    for entry_index = 1:numel(path_entries)
        current = path_entries{entry_index};
        if isempty(current)
            continue;
        end
        if strcmp(current, '.')
            current = pwd;
        end
        if ~isfolder(current)
            continue;
        end
        current = absolute_folder(current);

        while true
            [parent, folder_name] = fileparts(current);
            if strcmpi(folder_name, 'ppdata')
                [candidates, candidate_keys] = append_unique_folder( ...
                    candidates, candidate_keys, current);
            end

            child = fullfile(current, 'ppdata');
            if isfolder(child)
                [candidates, candidate_keys] = append_unique_folder( ...
                    candidates, candidate_keys, child);
            end

            if isempty(parent) || strcmp(parent, current)
                break;
            end
            current = parent;
        end
    end
end

function folder = absolute_folder(folder)
    [status, attributes] = fileattrib(folder);
    if status
        folder = attributes.Name;
    end
end

function [folders, keys] = append_unique_folder(folders, keys, folder)
    folder = absolute_folder(folder);
    key = folder;
    if ispc
        key = lower(key);
    end
    if ~any(strcmp(keys, key))
        folders{end + 1} = folder;
        keys{end + 1} = key;
    end
end

function filename = find_matching_pao_upf(root, symbol, funct, is_soc)
files = [dir(fullfile(root, '**', '*.UPF')); ...
    dir(fullfile(root, '**', '*.upf'))];
if isempty(files)
    error('compute_pao_projectability:NoUPFFiles', ...
        '赝势目录中没有 UPF 文件：%s', root);
end

functional_key = canonical_functional(funct);
best_score = -inf;
filename = '';
for ifile = 1:numel(files)
    [~, base, ext] = fileparts(files(ifile).name);
    if isempty(regexpi([base ext], ['^' regexptranslate('escape', symbol) ...
            '[._-]'], 'once'))
        continue;
    end

    candidate = fullfile(files(ifile).folder, files(ifile).name);
    text = fileread(candidate);
    if isempty(regexp(text, '<PP_CHI\.\d+[^>]*>[\s\S]*?\S[\s\S]*?</PP_CHI\.\d+>', ...
            'once'))
        continue;
    end

    lower_name = lower(candidate);
    score = 0;
    if contains(lower_name, functional_key)
        score = score + 100;
    end
    header_functional = regexp(text, ...
        'functional\s*=\s*"([^"]+)"', 'tokens', 'once');
    if ~isempty(header_functional) && ...
            functional_matches(header_functional{1}, functional_key)
        score = score + 50;
    end
    if contains(lower_name, 'default_fr')
        score = score + 5 * is_soc - 5 * ~is_soc;
    end
    if contains(lower(text(1:min(end, 8000))), 'is_paw="true"') || ...
            contains(lower(text(1:min(end, 8000))), 'is_ultrasoft="true"')
        score = score - 20;
    end

    if score > best_score
        best_score = score;
        filename = candidate;
    end
end

if isempty(filename)
    error('compute_pao_projectability:NoPAOInSelectedLibrary', ...
        ['The kssolvpptype-selected library "%s" contains no UPF with ', ...
         'nonempty PP_CHI.* data for element %s. PAO projectability ', ...
         'cannot be computed from beta projectors or default values.'], ...
        root, symbol);
end

if best_score < 50
    warning('compute_pao_projectability:FunctionalMismatchUsingSelectedPp', ...
        ['No PP_CHI UPF in the kssolvpptype-selected library "%s" matched ', ...
         'the requested/loaded functional "%s" for element %s. Using the ', ...
         'selected SCF pseudopotential file with PP_CHI data instead: %s'], ...
        root, funct, symbol, filename);
end

if isempty(filename)
    error('compute_pao_projectability:NoMatchingPAO', ...
        ['没有找到元素 %s、泛函 %s 对应且包含 PP_CHI.* 的 UPF。', ...
         '可投影性无法从 beta projector 代替计算。'], symbol, funct);
end
end

function key = canonical_functional(funct)
name = lower(strtrim(char(funct)));
if contains(name, 'pbe') || contains(name, 'hse')
    key = 'pbe';
elseif contains(name, 'pz') || contains(name, 'lda')
    key = 'pz';
elseif contains(name, 'blyp')
    key = 'blyp';
else
    key = name;
end
end

function tf = functional_matches(header, key)
header = lower(header);
switch key
    case 'pbe'
        tf = contains(header, 'pbe') || ...
            (contains(header, 'pbx') && contains(header, 'pbc'));
    case 'pz'
        tf = contains(header, 'pz') || contains(header, 'lda');
    otherwise
        tf = contains(header, key);
end
end

function [funct, source] = infer_pseudopotential_functional(crystal, crystal_funct)
%INFER_PSEUDOPOTENTIAL_FUNCTIONAL Prefer the UPF functional stored in ppvar.
    funct = strtrim(char(crystal_funct));
    source = 'crystal.funct';

    ppvar = read_crystal_member(crystal, 'ppvar', false);
    if isempty(ppvar)
        return;
    end

    pp_funct = [];
    try
        if isstruct(ppvar) && isfield(ppvar, 'funct')
            pp_funct = ppvar.funct;
        elseif isobject(ppvar) && isprop(ppvar, 'funct')
            pp_funct = ppvar.funct;
        end
    catch
        pp_funct = [];
    end

    if isempty(pp_funct)
        return;
    end
    if isstring(pp_funct)
        pp_funct = char(pp_funct);
    end
    if iscell(pp_funct)
        pp_funct = pp_funct{1};
    end
    if ischar(pp_funct) && ~isempty(strtrim(pp_funct))
        funct = strtrim(pp_funct);
        source = 'crystal.ppvar.funct';
    end
end

function [symbols, alist, xyzlist] = crystal_atom_data(crystal)
atoms = read_crystal_member(crystal, 'atoms', true);
alist = double(read_crystal_member(crystal, 'alist', true));
alist = alist(:);
xyzlist = double(read_crystal_member(crystal, 'xyzlist', true));

if size(xyzlist, 2) ~= 3 || numel(alist) ~= size(xyzlist, 1)
    error('compute_pao_projectability:BadAtoms', ...
        'crystal.xyzlist 与 crystal.alist 的尺寸不一致。');
end

symbols = cell(numel(alist), 1);
for iatom = 1:numel(alist)
    atom = atoms(alist(iatom));
    symbols{iatom} = strtrim(char(atom.symbol));
end
end

function value = read_crystal_member(crystal, name, required)
value = [];
try
    if isstruct(crystal) && isfield(crystal, name)
        value = crystal.(name);
    elseif isobject(crystal) && isprop(crystal, name)
        value = crystal.(name);
    end
catch
    value = [];
end
if isempty(value) && required
    error('compute_pao_projectability:MissingCrystalProperty', ...
        'crystal 缺少必需属性 %s。', name);
end
end
