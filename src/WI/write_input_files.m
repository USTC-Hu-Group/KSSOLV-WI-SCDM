function write_input_files(crystal, Eigvals, varargin)
%WRITE_INPUT_FILES Generate scdm.in and scdm.eig for the SCDM calculation.
%
% Supported calls:
%   write_input_files(crystal, Eigvals)
%   write_input_files(crystal, Eigvals, n_wannier, 'isolated')
%   write_input_files(crystal, Eigvals, 'erfc', mu, sigma)
%   write_input_files(crystal, Eigvals, n_wannier, 'erfc', mu, sigma)
%   write_input_files(crystal, Eigvals, [], 'gaussian', mu, sigma)
%   write_input_files(crystal, Eigvals, struct('num_wann', 8), ...)
%
% Eigvals is the KSSOLV eigenvalue matrix in Hartree. Its expected shape is
% [nbnd, nk]. Projectabilities are computed internally from PAOs stored in
% UPF pseudopotentials and Bloch wavefunctions stored in ./tmp/UNK* files.

    provenance = struct('status', 'started', 'stage', 'argument_parsing', ...
        'created_at', char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), ...
        'working_directory', pwd, 'arguments', {varargin});
    try

    % Preserve old positional calls while allowing n_wannier to be omitted.
    [n_wannier, focc, mu, sigma, has_manual_mu_sigma] = ...
        parse_write_input_args(varargin);
    focc = normalize_entanglement_mode(focc);
    provenance.mode = focc;
    provenance.manual_parameters = has_manual_mu_sigma;
    provenance.manual_parameters_sanitized = false;
    provenance.sanitized_fields = {};
    provenance.parameter_source = 'automatic_fit';

    if strcmp(focc, 'isolated')
        provenance.parameter_source = 'not_applicable';
        % Isolated bands use uniform weights; mu and sigma are irrelevant.
        mu = [];
        sigma = [];
        manual_mu = [];
        manual_sigma = [];
    else
        % Manual values are used only when both mu and sigma are supplied.
        if has_manual_mu_sigma
            provenance.requested_mu = mu;
            provenance.requested_sigma = sigma;
            provenance.parameter_source = 'manual';
        end
        [manual_mu, manual_sigma, sanitized_fields] = sanitize_manual_parameters(mu, sigma);
        if has_manual_mu_sigma
            provenance.sanitized_fields = sanitized_fields;
            provenance.manual_parameters_sanitized = ~isempty(sanitized_fields);
            if provenance.manual_parameters_sanitized
                provenance.parameter_source = 'manual_sanitized';
            end
        end
    end

    output = 'scdm.in';
    nk = crystal.nkpts;
    crystal_nbnd = crystal.nbnd;
    kpts = crystal.kpts;

    validateattributes(nk, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(crystal_nbnd, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});
    if size(kpts, 1) ~= nk || size(kpts, 2) < 3
        error('write_input_files:BadKpoints', ...
            'crystal.kpts must contain nk rows and at least three columns.');
    end
    kpts = kpts(:, 1:3);

    % nbnd is determined from the actual eigenvalue data, not crystal.nbnd.
    % This permits, for example, 12 input bands with 8 target Wannier
    % orbitals even if crystal.nbnd has already been changed to 8.
    nbnd = infer_band_count(Eigvals, nk);
    if crystal_nbnd ~= nbnd
        warning('write_input_files:CrystalBandCountMismatch', ...
            ['crystal.nbnd=%d，但 Eigvals 包含 %d 条能带。', ...
             'scdm.eig 将保留 Eigvals 中的全部 %d 条能带。'], ...
            crystal_nbnd, nbnd, nbnd);
    end

    % KSSOLV eigenvalues are stored in Hartree; all SCDM energy parameters
    % and the generated scdm.eig file use eV.
    eig_hartree = normalize_band_k_matrix(Eigvals, nbnd, nk, 'Eigvals');
    eigs = eig_hartree * 27.211386246;
    provenance.energy_eV = eigs;
    provenance.energy_reference = 'absolute KSSOLV eigenvalues; Hartree times 27.211386246';
    provenance.kpoints_cartesian = kpts;
    provenance.cell_bohr = double(crystal.supercell);
    provenance.atoms_bohr = double(crystal.xyzlist);
    provenance.stage = 'pao_construction';

    [has_manual_num_wann, manual_num_wann] = ...
        parse_num_wann_override(n_wannier);

    % PAO 构造按实际需求执行：
    % 1. 未手动指定 mu/sigma 的 erfc/gaussian 模式：计算完整可投影性；
    % 2. 仅 num_wann 需要自动判定：只读取 UPF 统计 N_PAO，不读取 UNK；
    % 3. num_wann 和 mu/sigma 均已手动指定：完全跳过 PAO/UNK 处理。
    proj_matrix = [];
    pao_info = struct();
    needs_auto_fit = ~strcmp(focc, 'isolated') && ~has_manual_mu_sigma;
    needs_pao_count = ~has_manual_num_wann;
    provenance.manual_parameters = has_manual_mu_sigma;
    provenance.manual_num_wann = has_manual_num_wann;
    if needs_auto_fit || needs_pao_count
        count_only = ~needs_auto_fit;
        [proj_matrix, pao_info] = compute_pao_projectability( ...
            crystal, nbnd, nk, count_only);
    end
    provenance.pao_info = pao_info;
    provenance.projectability = proj_matrix;
    provenance.stage = 'target_dimension';

    %% === START OF AUTO num_wann DETERMINATION ===
    % Vitale et al. (2020) 高通量协议令目标 Wannier 轨道数 J 等于
    % 构造投影时采用的 PAO 总数 N_PAO。这里采用“手动优先、自动保底”：
    %
    % 1. varargin 中若显式给出正整数，或 options.num_wann/options.n_wannier
    %    是正整数，则严格使用该手动值，不被自动结果覆盖。
    % 2. 若第三个输入为空、为零或为负数，则从 PAO 构造模块返回的
    %    pao_info.total_count 中获取 N_PAO。SOC 情况包含两个旋量分量。
    % 3. PAO 构造失败、计数无效或 N_PAO 超过当前 nbnd 时直接报错。
    %    这避免在批量测试中静默写入不可信的 num_wann。
    if has_manual_num_wann
        n_wannier = manual_num_wann;
        if n_wannier > nbnd
            error('write_input_files:TooManyWannierFunctions', ...
                '手动 num_wann=%d 不能超过输入能带数 nbnd=%d。', ...
                n_wannier, nbnd);
        end
        fprintf('num_wann 使用手动重载值：%d\n', n_wannier);
    else
        if ~isfield(pao_info, 'total_count')
            error('write_input_files:MissingPAOCount', ...
                ['PAO projectability did not return total_count; ', ...
                 'cannot determine num_wann automatically.']);
        end
        auto_num_wann = pao_info.total_count;
        validateattributes(auto_num_wann, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'integer', 'positive'});
        if auto_num_wann > nbnd
            error('write_input_files:PAOCountExceedsBands', ...
                ['自动得到 N_PAO=%d，但当前仅计算了 nbnd=%d 条能带；', ...
                 '无法构造多于能带数的 Wannier 轨道。'], ...
                auto_num_wann, nbnd);
        end
        n_wannier = double(auto_num_wann);
        fprintf(['num_wann 自动设为 N_PAO：%d（空间 PAO=%d，', ...
            '旋量倍数=%d）\n'], n_wannier, ...
            pao_info.spatial_count, pao_info.spin_factor);
    end
    %% === END OF AUTO num_wann DETERMINATION ===
    provenance.num_wann = n_wannier;
    provenance.stage = 'parameter_selection';

    %% === START OF AUTO-SCDM PARAMETERS (Vitale et al. 2020) ===
    % 可投影性 p_nk 衡量 Bloch 态在局域 PAO 子空间中的权重。局域性良好的
    % 低能态通常满足 p_nk ~ 1，而高能自由电子态通常逐渐衰减至 p_nk ~ 0。
    %
    % 将所有能带和 k 点展平后，用互补误差函数
    %   p(E) = 1/2 erfc((E - mu_fit) / sigma_fit)
    % 拟合这一能量依赖关系。mu_fit 描述可投影性下降中心，sigma_fit 描述
    % 过渡区的能量宽度。根据 Vitale et al. (2020) 的参数选取方式，最终
    % 写入 SCDM 的参数为 mu = mu_fit - 3*sigma_fit，sigma = sigma_fit。
    %
    % compute_pao_projectability 从 MATLAB 路径中定位 ppdata，通过
    % kssolvpptype 锁定子库，再按 crystal.funct 选择 UPF/PP_CHI.*
    % 构造 PAO，并从 ./tmp/UNK* 读取 Bloch 波函数。PAO 在每个 k 点通过
    % SVD 正交化后形成投影算符，因此得到的 p_nk 严格位于 [0,1]。
    %
    % isolated 模式使用恒等过滤，不计算 mu/sigma。其余模式若已手动指定
    % mu/sigma，则严格使用手动值；未手动指定时才执行自动拟合。
    if strcmp(focc, 'isolated')
        fprintf('SCDM 模式为 isolated：跳过 mu 与 sigma 的计算。\n');
    elseif has_manual_mu_sigma
        mu = manual_mu;
        sigma = manual_sigma;
        fprintf(['使用手动 SCDM 参数：mu = %.8f eV, ', ...
            'sigma = %.8f eV；跳过可投影性拟合。\n'], mu, sigma);
    else
        proj_matrix = normalize_band_k_matrix( ...
            proj_matrix, nbnd, nk, 'computed projectability');

        energies_fit = eigs(:);
        projectability_fit = double(proj_matrix(:));
        valid = isfinite(energies_fit) & isfinite(projectability_fit);
        provenance.valid_fit_mask = reshape(valid, nbnd, nk);
        energies_fit = energies_fit(valid);
        projectability_fit = projectability_fit(valid);

        if numel(energies_fit) < 4
            error('write_input_files:TooFewFitData', ...
                '有效的能量-可投影性数据点少于 4 个。');
        end
        if any(projectability_fit < -1e-8 | projectability_fit > 1 + 1e-8)
            error('write_input_files:BadProjectabilityRange', ...
                '自动计算的可投影性中存在明显超出 [0,1] 的数值。');
        end

        % 容忍浮点计算造成的微小越界，但拒绝真正不物理的脏数据。
        projectability_fit = min(max(projectability_fit, 0), 1);
        if max(projectability_fit) - min(projectability_fit) < 0.05
            error('write_input_files:FlatProjectability', ...
                '自动计算的可投影性变化范围过小，无法确定阶跃衰减位置。');
        end

        % mu_fit 初猜取 p_nk 最接近 0.5 的态对应的能量；这是阶跃中心
        % 最直接且物理上稳定的估计。sigma_fit 初猜固定为 2.0 eV。
        [~, half_index] = min(abs(projectability_fit - 0.5));
        mu_fit_guess = energies_fit(half_index);
        sigma_fit_guess = 2.0;

        energy_min = min(energies_fit);
        energy_max = max(energies_fit);
        energy_span = max(energy_max - energy_min, 1.0);

        % sigma_fit 必须具有正的物理展宽意义，因此严格限制为 >0.1 eV。
        sigma_lower_bound = 0.100001;
        lower_bound = [energy_min - energy_span, sigma_lower_bound];
        upper_bound = [energy_max + energy_span, max(2 * energy_span, 20.0)];
        initial_guess = [mu_fit_guess, sigma_fit_guess];
        initial_guess = min(max(initial_guess, lower_bound), upper_bound);
        provenance.fit_initial_guess = initial_guess;
        provenance.fit_lower_bound = lower_bound;
        provenance.fit_upper_bound = upper_bound;

        erfc_model = @(x, energy) ...
            0.5 * erfc((energy - x(1)) ./ x(2));
        fit_options = optimoptions('lsqcurvefit', ...
            'Display', 'off', ...
            'FunctionTolerance', 1e-12, ...
            'StepTolerance', 1e-12, ...
            'OptimalityTolerance', 1e-10, ...
            'MaxFunctionEvaluations', 5000);

        [fit_result, fit_resnorm, fit_residual, exitflag, fit_output] = lsqcurvefit( ...
            erfc_model, initial_guess, energies_fit, ...
            projectability_fit, lower_bound, upper_bound, fit_options);
        provenance.fit_result = fit_result;
        provenance.fit_resnorm = fit_resnorm;
        provenance.fit_residual = fit_residual;
        provenance.fit_exitflag = exitflag;
        provenance.fit_output = fit_output;
        provenance.fit_rmse = sqrt(mean(fit_residual.^2));

        mu_fit = fit_result(1);
        sigma_fit = fit_result(2);
        fitted_mu = mu_fit - 3 * sigma_fit;
        fitted_sigma = sigma_fit;

        if exitflag <= 0 || ~isreal(fitted_mu) || ~isreal(fitted_sigma) || ...
                ~isfinite(fitted_mu) || ~isfinite(fitted_sigma) || ...
                fitted_sigma <= 0.1
            error('write_input_files:UnreliableFit', ...
                '非线性拟合未可靠收敛，或返回了非物理参数。');
        end

        % 在写盘前直接覆写原有局部变量。
        mu = fitted_mu;
        sigma = fitted_sigma;
        fprintf(['Auto-SCDM 参数拟合完成: mu_fit = %.8f eV, ', ...
            'sigma_fit = %.8f eV, mu = %.8f eV, sigma = %.8f eV\n'], ...
            mu_fit, sigma_fit, mu, sigma);
    end
    %% === END OF AUTO-SCDM PARAMETERS (Vitale et al. 2020) ===
    provenance.mu_eV = mu;
    provenance.sigma_eV = sigma;
    provenance.stage = 'file_output';

    % Crystal.finalize stores k-points in Cartesian reciprocal coordinates.
    % Convert them back to coefficients of the reciprocal lattice basis.
    [k_frac, k_grid_dims] = crystal_kpoints_to_fractional(crystal, kpts, nk);

    fprintf('正在生成 %s ...\n', output);
    fid = fopen(output, 'w');
    if fid == -1
        error('write_input_files:CannotOpenInput', ...
            '无法创建文件 %s。', output);
    end
    input_cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '%d %d %d\n', nbnd, nk, n_wannier);
    if strcmp(focc, 'isolated')
        fprintf(fid, 'scdm.eig %s\n', focc);
    else
        fprintf(fid, 'scdm.eig %s %.12f %.12f\n', focc, mu, sigma);
    end
    for ik = 1:nk
        fprintf(fid, '%.10f %.10f %.10f\n', ...
            k_frac(ik, 1), k_frac(ik, 2), k_frac(ik, 3));
    end
    write_cube_metadata(fid, crystal, k_grid_dims, pao_info);
    clear input_cleanup

    fid = fopen('scdm.eig', 'w');
    if fid == -1
        error('write_input_files:CannotOpenEigenvalues', ...
            '无法创建文件 scdm.eig。');
    end
    eig_cleanup = onCleanup(@() fclose(fid));

    % 外层 k 点、内层能带，与 Wannier90/SCDM 的 Fortran 读取顺序一致。
    for ik = 1:nk
        for ibnd = 1:nbnd
            fprintf(fid, '%d %d %.12f\n', ibnd, ik, eigs(ibnd, ik));
        end
    end
    clear eig_cleanup

    fprintf('成功生成 SCDM 输入文件（k 点数: %d，能带数: %d）。\n', nk, nbnd);
    provenance.status = 'success';
    provenance.stage = 'complete';
    catch failure
        provenance.status = 'failed';
        provenance.error_identifier = failure.identifier;
        provenance.error_message = failure.message;
        save_scdm_provenance(provenance);
        rethrow(failure);
    end
    save_scdm_provenance(provenance);
end

function write_cube_metadata(fid, crystal, k_grid_dims, pao_info)
%WRITE_CUBE_METADATA Append structure data needed for Gaussian Cube output.
    if nargin < 4
        pao_info = struct();
    end

    supercell = double(crystal.supercell);
    xyzlist = double(crystal.xyzlist);
    alist = double(crystal.alist(:));
    atoms = crystal.atoms;

    if ~isequal(size(supercell), [3, 3])
        error('write_input_files:BadSupercell', ...
            'crystal.supercell must be a 3-by-3 matrix in Bohr.');
    end
    if size(xyzlist, 2) ~= 3 || size(xyzlist, 1) ~= numel(alist)
        error('write_input_files:BadAtomCoordinates', ...
            'crystal.xyzlist must contain one Cartesian coordinate per atom.');
    end

    fprintf(fid, 'BEGIN_CUBE_METADATA\n');
    fprintf(fid, 'KGRID %d %d %d\n', k_grid_dims);
    fprintf(fid, 'CELL_BOHR\n');
    for row = 1:3
        fprintf(fid, '%.16e %.16e %.16e\n', supercell(row, :));
    end

    fprintf(fid, 'ATOMS_BOHR %d\n', numel(alist));
    for atom_index = 1:numel(alist)
        atom_type = alist(atom_index);
        atom = atoms(atom_type);
        symbol = strtrim(char(atom.symbol));
        atomic_number = read_finite_atom_field(atom, 'anum');
        effective_charge = read_finite_atom_field(atom, 'venum');
        if ~isfinite(effective_charge)
            effective_charge = read_z_valence_from_pao_info(symbol, pao_info);
            if isfinite(effective_charge)
                fprintf(['Cube metadata effective charge for %s read ', ...
                    'from UPF z_valence: %.8g\n'], ...
                    symbol, effective_charge);
            end
        end
        if ~isfinite(atomic_number) || ~isfinite(effective_charge)
            error('write_input_files:MissingCubeAtomMetadata', ...
                ['Cannot write Cube metadata for element %s: atomic number ', ...
                 'or effective charge is not finite. Check Atom data and ', ...
                 'the selected UPF z_valence.'], symbol);
        end
        fprintf(fid, '%d %.10f %.16e %.16e %.16e\n', ...
            atomic_number, effective_charge, ...
            xyzlist(atom_index, :));
    end
    fprintf(fid, 'END_CUBE_METADATA\n');
end

function value = read_finite_atom_field(atom, field_name)
    value = NaN;
    try
        candidate = double(atom.(field_name));
        if isscalar(candidate) && isfinite(candidate)
            value = candidate;
        end
    catch
        value = NaN;
    end
end

function z_valence = read_z_valence_from_pao_info(symbol, pao_info)
    z_valence = NaN;
    if ~isstruct(pao_info) || ~isfield(pao_info, 'pp_files')
        return;
    end

    pp_files = pao_info.pp_files;
    if isstring(pp_files)
        pp_files = cellstr(pp_files);
    end
    if ~iscell(pp_files)
        return;
    end

    for ifile = 1:numel(pp_files)
        filename = pp_files{ifile};
        if isempty(filename) || ~(ischar(filename) || isstring(filename))
            continue;
        end
        filename = char(filename);
        [~, base, ext] = fileparts(filename);
        if isempty(regexpi([base ext], ...
                ['^' regexptranslate('escape', symbol) '[._-]'], 'once'))
            continue;
        end
        z_valence = read_upf_z_valence(filename);
        if isfinite(z_valence)
            return;
        end
    end
end

function z_valence = read_upf_z_valence(filename)
    z_valence = NaN;
    if ~isfile(filename)
        return;
    end
    text = fileread(filename);
    token = regexpi(text, 'z_valence\s*=\s*"([^"]+)"', ...
        'tokens', 'once');
    if isempty(token)
        token = regexpi(text, 'z_valence\s*=\s*''([^'']+)''', ...
            'tokens', 'once');
    end
    if isempty(token)
        token = regexpi(text, 'z_valence\s*=\s*([-+0-9.eEdD]+)', ...
            'tokens', 'once');
    end
    if isempty(token)
        return;
    end
    z_valence = str2double(strrep(strtrim(token{1}), 'D', 'E'));
end

function [k_frac, k_grid_dims] = crystal_kpoints_to_fractional(crystal, kpts, nk)
%CRYSTAL_KPOINTS_TO_FRACTIONAL Normalize KSSOLV k-points for SCDM.

    supercell = double(crystal.supercell);
    reciprocal_cell = 2 * pi * inv(supercell).';
    k_frac = double(kpts) / reciprocal_cell;
    k_frac(abs(k_frac) < 1e-12) = 0;

    k_grid_dims = [];
    try
        k_grid_dims = double(crystal.nkxyz);
    catch
    end
    k_grid_dims = k_grid_dims(:).';

    if numel(k_grid_dims) ~= 3 || any(~isfinite(k_grid_dims)) || ...
            any(k_grid_dims < 1) || any(k_grid_dims ~= round(k_grid_dims)) || ...
            prod(k_grid_dims) ~= nk
        k_grid_dims = infer_fractional_grid_dims(k_frac, nk);
    end

    validate_fractional_kgrid(k_frac, k_grid_dims, nk);
    fprintf(['Converted KSSOLV Cartesian k-points to fractional ', ...
        'coordinates; KGRID = %d %d %d.\n'], ...
        k_grid_dims(1), k_grid_dims(2), k_grid_dims(3));
end

function dims = infer_fractional_grid_dims(k_frac, nk)
    dims = zeros(1, 3);
    rounded = round(k_frac * 1e9) / 1e9;
    for axis = 1:3
        dims(axis) = numel(unique(rounded(:, axis)));
    end
    if prod(dims) ~= nk
        error('write_input_files:CannotInferKGrid', ...
            ['Cannot infer a complete uniform k-grid from fractional ', ...
             'coordinates: inferred %d x %d x %d, nk=%d.'], ...
            dims(1), dims(2), dims(3), nk);
    end
end

function validate_fractional_kgrid(k_frac, dims, nk)
    scaled = k_frac .* dims;
    offset = scaled(1, :) - round(scaled(1, :));
    integer_coords = round(scaled - offset);
    residual = max(abs(scaled - offset - integer_coords), [], 'all');
    indices = mod(integer_coords, dims) + 1;
    linear_indices = sub2ind(dims, indices(:, 1), indices(:, 2), indices(:, 3));

    if residual > 1e-6 || numel(unique(linear_indices)) ~= nk
        error('write_input_files:InvalidUniformKGrid', ...
            ['The k-points do not form a complete uniform %d x %d x %d ', ...
             'grid; maximum integer-grid residual is %.3e.'], ...
            dims(1), dims(2), dims(3), residual);
    end
end

function nbnd = infer_band_count(values, nk)
%INFER_BAND_COUNT Determine the number of input bands from eigenvalues.
    validateattributes(values, {'numeric'}, {'2d', 'nonempty', 'real'}, ...
        'write_input_files', 'Eigvals');

    if isvector(values)
        if mod(numel(values), nk) ~= 0
            error('write_input_files:BadEigenvalueCount', ...
                'Eigvals 长度 %d 不能被 k 点数 nk=%d 整除。', ...
                numel(values), nk);
        end
        nbnd = numel(values) / nk;
    elseif size(values, 2) == nk
        nbnd = size(values, 1);
    elseif size(values, 1) == nk
        nbnd = size(values, 2);
    else
        error('write_input_files:BadEigenvalueSize', ...
            ['Eigvals 必须为 [nbnd,nk]、[nk,nbnd]，或长度为 ', ...
             'nbnd*nk 的向量；当前尺寸为 [%d,%d]，nk=%d。'], ...
            size(values, 1), size(values, 2), nk);
    end

    validateattributes(nbnd, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});
end

function [n_wannier, focc, mu, sigma, has_manual_mu_sigma] = ...
        parse_write_input_args(args)
%PARSE_WRITE_INPUT_ARGS Parse optional n_wannier/focc/mu/sigma inputs.
    n_wannier = [];
    focc = 'isolated';
    mu = 0.0;
    sigma = 0.1;
    has_manual_mu_sigma = false;

    if isempty(args)
        return;
    end

    if numel(args) == 1 && isstruct(args{1})
        options = args{1};
        n_wannier = read_option_field(options, {'num_wann', 'n_wannier'}, []);
        focc = read_option_field(options, {'focc', 'entanglement'}, focc);
        mu = read_option_field(options, {'mu'}, mu);
        sigma = read_option_field(options, {'sigma'}, sigma);
        has_manual_mu_sigma = isfield(options, 'mu') && ...
            isfield(options, 'sigma');
        return;
    end

    if numel(args) == 1
        if ischar(args{1}) || (isstring(args{1}) && isscalar(args{1}))
            focc = args{1};
        else
            n_wannier = args{1};
        end
        return;
    end

    % Isolated mode needs no mu or sigma.
    %   write_input_files(crystal, Eigvals, n_wannier, 'isolated')
    if numel(args) == 2
        n_wannier = args{1};
        focc = args{2};
        return;
    end

    % New shorthand without n_wannier:
    %   write_input_files(crystal, Eigvals, 'erfc', mu, sigma)
    if numel(args) == 3
        focc = args{1};
        mu = args{2};
        sigma = args{3};
        has_manual_mu_sigma = true;
        return;
    end

    % Backward-compatible positional form:
    %   write_input_files(crystal, Eigvals, n_wannier, 'erfc', mu, sigma)
    if numel(args) == 4
        n_wannier = args{1};
        focc = args{2};
        mu = args{3};
        sigma = args{4};
        has_manual_mu_sigma = true;
        return;
    end

    error('write_input_files:BadOptionalInputs', ...
        ['可选参数应为 (), (mode), (n_wannier, isolated), ', ...
         '(mode, mu, sigma)，或 (n_wannier, mode, mu, sigma)。']);
end

function value = read_option_field(options, names, default_value)
%READ_OPTION_FIELD Read the first existing field from a struct.
    value = default_value;
    for iname = 1:numel(names)
        if isfield(options, names{iname})
            value = options.(names{iname});
            return;
        end
    end
end

function mode = normalize_entanglement_mode(value)
%NORMALIZE_ENTANGLEMENT_MODE Validate the named SCDM filtering mode.
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('write_input_files:NumericFoccNotSupported', ...
            ['SCDM 过滤模式必须使用字符输入：''isolated''、''erfc'' ', ...
             '或 ''gaussian''；不再接受数字 focc。']);
    end
    mode = lower(strtrim(char(value)));
    valid_modes = {'isolated', 'erfc', 'gaussian'};
    if ~any(strcmp(mode, valid_modes))
        error('write_input_files:BadEntanglementMode', ...
            '未知 SCDM 过滤模式 %s。可选值为 isolated、erfc、gaussian。', ...
            mode);
    end
end

function [has_override, num_wann] = parse_num_wann_override(value)
%PARSE_NUM_WANN_OVERRIDE Extract a valid explicit num_wann setting.
    has_override = false;
    num_wann = [];

    if isstruct(value)
        if isfield(value, 'num_wann')
            value = value.num_wann;
        elseif isfield(value, 'n_wannier')
            value = value.n_wannier;
        else
            value = [];
        end
    end

    if isempty(value) || ...
            (isnumeric(value) && isscalar(value) && isreal(value) && ...
             isfinite(value) && value <= 0)
        return;
    end

    if isnumeric(value) && isscalar(value) && isreal(value) && ...
            isfinite(value) && value > 0 && value == round(value)
        has_override = true;
        num_wann = double(value);
        return;
    end

    warning('write_input_files:InvalidManualNumWann', ...
        '手动 num_wann 不是有效正整数，将忽略该值并尝试自动判定。');
end

function matrix = normalize_band_k_matrix(values, nbnd, nk, variable_name)
%NORMALIZE_BAND_K_MATRIX Return data in the canonical [nbnd, nk] layout.
    validateattributes(values, {'numeric'}, {'2d', 'nonempty', 'real'}, ...
        'write_input_files', variable_name);

    if isequal(size(values), [nbnd, nk])
        matrix = double(values);
    elseif isequal(size(values), [nk, nbnd]) && nbnd ~= nk
        matrix = double(values.');
    elseif isvector(values) && numel(values) == nbnd * nk
        matrix = reshape(double(values), nbnd, nk);
    else
        error('write_input_files:BadMatrixSize', ...
            '%s 的尺寸必须为 [%d,%d]，或包含相同数量的线性数据。', ...
            variable_name, nbnd, nk);
    end
end

function [mu, sigma, sanitized_fields] = sanitize_manual_parameters(mu, sigma)
%SANITIZE_MANUAL_PARAMETERS Preserve legacy replacements, not physical guarantees.
    sanitized_fields = {};
    if ~(isnumeric(mu) && isscalar(mu) && isreal(mu) && isfinite(mu))
        warning('write_input_files:BadManualMu', ...
            '手动 mu 无效，使用兼容性替代值 0.0 eV；该值不保证物理适用性。');
        mu = 0.0;
        sanitized_fields{end+1} = 'mu';
    else
        mu = double(mu);
    end

    if ~(isnumeric(sigma) && isscalar(sigma) && isreal(sigma) && ...
            isfinite(sigma) && sigma > 0)
        warning('write_input_files:BadManualSigma', ...
            '手动 sigma 无效，使用兼容性替代值 2.0 eV；该值不保证物理适用性。');
        sigma = 2.0;
        sanitized_fields{end+1} = 'sigma';
    else
        sigma = double(sigma);
    end
end
