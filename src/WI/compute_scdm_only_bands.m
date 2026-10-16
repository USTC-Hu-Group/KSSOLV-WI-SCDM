function [interpolated_bands, H_R, transform_info] = compute_scdm_only_bands( ...
    Amat_in, eigenvalues_in, kpts_in, k_path_fine, efermi, varargin)
% COMPUTE_SCDM_ONLY_BANDS - 计算并绘制 "SCDM-only" 插值能带
% 输入参数：
%   Amat_in            - 投影矩阵文件，名为"xx.amn"
%   eigenvalues.in     - NDFT 本征值文件
%   kpts.in            -k点列表
%   k_path_fine     - [n_fine_pts x 3] 目标插值路径的分数坐标
%
% 输出参数：
%   interpolated_bands - [n_fine_pts x n_wannier] 插值后的能带 (单位: eV)
%   H_R                - 实空间哈密顿量 (单位: eV)
%                        包含网格原点相位，使用居中的整数 R 索引。
%   transform_info     - 变换参数及 kgrid_dims / kgrid_shift_fractional。
%% 1. 常数定义与参数检查

options = parse_interpolation_options(varargin{:});
use_wi_scdm_f = strcmp(options.Method, 'wi-scdm-f');

%Amat_in 处理
famn = fopen(Amat_in,"r");
line = fgetl(famn);
disp(line);
nband = fscanf(famn,'%d',1);
nk = fscanf(famn,'%d',1);
n_wannier = fscanf(famn,'%d',1);
data = fscanf(famn, '%d %d %d %f %f', [5, Inf])';
fclose(famn);
nband_idx = data(:,1);
n_wannier_idx = data(:,2);
kpt_idx = data(:,3);
real_part = data(:,4);
imag_part = data(:,5);
Amat = zeros(nband,n_wannier,nk);
lin_idx = sub2ind([nband, n_wannier, nk], nband_idx, n_wannier_idx, kpt_idx);
Amat(lin_idx) = real_part + 1i * imag_part;

%eig处理
feig = fopen(eigenvalues_in,"r");
data_1 = fscanf(feig,'%d %d %f',[3,Inf])';
fclose(feig);
eigenvalue = data_1(:,3);
E_bands_eV = reshape(eigenvalue,nband,nk);

if use_wi_scdm_f
    if isempty(options.TransformParams)
        transform_info = wi_scdm_f_params(E_bands_eV, ...
            'TransformN', options.TransformN, ...
            'TransformA', options.TransformA, ...
            'TransformEpsilon', options.TransformEpsilon, ...
            'TransformTopBand', options.TransformTopBand, ...
            'InverseSamples', options.InverseSamples);
    else
        transform_info = options.TransformParams;
        transform_info.method = 'wi-scdm-f';
        if ~isfield(transform_info, 'top_band') || isempty(transform_info.top_band)
            transform_info.top_band = NaN;
        end
        if ~isfield(transform_info, 'inverse_samples') || ...
                isempty(transform_info.inverse_samples)
            transform_info.inverse_samples = options.InverseSamples;
        end
    end
    E_for_hamiltonian_eV = wi_scdm_f_transform(E_bands_eV, transform_info);
else
    transform_info = struct('method', options.Method, ...
        'h_r_is_transformed', false);
    E_for_hamiltonian_eV = E_bands_eV;
end


%k点提取
[~, k_grid_dims, grid_linear_indices, grid_shift] = ...
    read_scdm_uniform_kgrid(kpts_in, nk);
transform_info.kgrid_dims = k_grid_dims;
transform_info.kgrid_shift_fractional = grid_shift;
transform_info.fourier_convention = 'centered-R-with-mesh-origin-phase-v1';

[n_fine, ~] = size(k_path_fine);

fprintf('=== 开始计算 SCDM 插值能带 ===\n');
fprintf('  - Wannier 函数数量: %d\n', n_wannier);
if use_wi_scdm_f
    transform_info.h_r_is_transformed = true;
    fprintf(['  - 插值方法: WI-SCDM-f, n = %.3g, a = %.8f eV, ', ...
        'epsilon = %.8f eV, top band = %d\n'], ...
        transform_info.n, transform_info.a, transform_info.epsilon, ...
        transform_info.top_band);
    fprintf(['  - 注意: 返回的 H_R 是 f(H) 的实空间矩阵，', ...
        '不是未变换的物理哈密顿量。\n']);
else
    fprintf('  - 插值方法: WI-SCDM/SCDM-only\n');
end

%% 2. 构建 Wannier 规范哈密顿量 H(k)
% H_wann(k) = U' * H_dft(k) * U
% 其中 U 由 SCDM 投影矩阵 A 经 Löwdin 正交化得到
H_wann_k = zeros(k_grid_dims(1), k_grid_dims(2), k_grid_dims(3), n_wannier, n_wannier);
fprintf('  正在执行 Löwdin 正交化和规范变换...\n');
for ik = 1:nk
    % 1. 获取投影矩阵 A [nbands x n_wannier]
    A_k = Amat(:, :, ik);
    
    % 2. Löwdin 正交化: U = A * (A'A)^(-1/2)
    S = A_k' * A_k;
    [V, D] = eig(S);
    
    % 防止数值噪音导致的除零
    d_diag = diag(D);
    d_diag(d_diag < 1e-10) = 1e-10;
    
    S_inv_sqrt = V * diag(1./sqrt(d_diag)) * V';
    U_k = A_k * S_inv_sqrt;
    
    % 3. 构建哈密顿量 (使用 eV 单位的能量)
    % E_bands_eV 是 [nbands x nkpts] 矩阵，直接取第 ik 列
    E_k_diag = diag(E_for_hamiltonian_eV(:, ik));
    
    H_local = U_k' * E_k_diag * U_k;
    
    % 4. 映射到 FFT 网格
    [idx1, idx2, idx3] = ind2sub(k_grid_dims, grid_linear_indices(ik));
    idx = [idx1, idx2, idx3];
    
    H_wann_k(idx(1), idx(2), idx(3), :, :) = H_local;
end
%% 3. 傅里叶变换到实空间 H(R)
[R1, R2, R3] = ndgrid(0:k_grid_dims(1)-1, 0:k_grid_dims(2)-1, 0:k_grid_dims(3)-1);
% Center the integer FFT indices (not a metric-based Wigner-Seitz cell).
R1(R1 >= k_grid_dims(1)/2) = R1(R1 >= k_grid_dims(1)/2) - k_grid_dims(1);
R2(R2 >= k_grid_dims(2)/2) = R2(R2 >= k_grid_dims(2)/2) - k_grid_dims(2);
R3(R3 >= k_grid_dims(3)/2) = R3(R3 >= k_grid_dims(3)/2) - k_grid_dims(3);
% For k_m = m./N + delta, FFT returns C_R relative to the mesh origin.
% The physical Fourier coefficients are H_R = exp(-2*pi*i*delta.R) C_R.
origin_phase = exp(-2i*pi*(R1.*grid_shift(1) + ...
    R2.*grid_shift(2) + R3.*grid_shift(3)));
fprintf('  执行 FFT 变换到实空间 H(R)...\n');
H_R = zeros(size(H_wann_k));
H_fft_flat = zeros(nk, n_wannier, n_wannier);
for i = 1:n_wannier
    for j = 1:n_wannier
        coefficients = fftn(H_wann_k(:, :, :, i, j)) / nk;
        H_R(:, :, :, i, j) = coefficients .* origin_phase;
        H_fft_flat(:, i, j) = coefficients(:);
    end
end
%% 4. 傅里叶插值回 k' 空间
interpolated_bands = zeros(n_fine, n_wannier);
fprintf('  正在插值到 %d 个精细 k 点...\n', n_fine);
for ik = 1:n_fine
    % C_R is evaluated at q = k' - delta, including all Nyquist terms.
    k_prime = k_path_fine(ik, :) - grid_shift;
    
    % 计算相位
    % H(k') = sum_R C(R) * exp(2*pi*i*(k'-delta).R).
    % 对居中的整数 FFT 网格进行边界对称化：
    % 对于偶数网格 N，边界点 -N/2 实际上对应 aliased 的 +N/2。
    % 为保证插值平滑性和 Hermiticity，需将边界项权重设为 0.5，并在 +N/2 处添加镜像项。
    % 简而言之，对于边界维度 d，将 e^(i k_d R_d) 替换为 cos(k_d * N_d * pi)。
    
    % 初始化相位因子
    % phase_factor = exp(2i * pi * (R1*k1 + R2*k2 + R3*k3))
    phase_exponent = R1(:)*k_prime(1) + R2(:)*k_prime(2) + R3(:)*k_prime(3);
    phase_factor = exp(2i * pi * phase_exponent);
    
    % 应用 FFT 网格边界修正 (Symmetrization)
    % 针对每个维度检查偶数网格的边界点
    
    % Dimension 1
    if mod(k_grid_dims(1), 2) == 0
        boundary_idx_1 = (R1(:) == -k_grid_dims(1)/2);
        % 将 e^(-i pi k N) 替换为 cos(pi k N)
        % 原相位因子中该维度的贡献是 exp(2pi * i * k * (-N/2)) = exp(-i * pi * k * N)
        % 我们将其修正为 0.5 * (exp(...) + exp(+...)) = cos(pi * k * N)
        % 实现方式：乘以 cos/exp
        
        correction_factor = cos(pi * k_prime(1) * k_grid_dims(1)) ./ ...
                            exp(2i * pi * k_prime(1) * (-k_grid_dims(1)/2));
        
        % 避免除零（理论上 exp 不会为0，但为了数值安全）
        phase_factor(boundary_idx_1) = phase_factor(boundary_idx_1) .* correction_factor;
    end
    
    % Dimension 2
    if mod(k_grid_dims(2), 2) == 0
        boundary_idx_2 = (R2(:) == -k_grid_dims(2)/2);
        correction_factor = cos(pi * k_prime(2) * k_grid_dims(2)) ./ ...
                            exp(2i * pi * k_prime(2) * (-k_grid_dims(2)/2));
        phase_factor(boundary_idx_2) = phase_factor(boundary_idx_2) .* correction_factor;
    end
    
    % Dimension 3
    if mod(k_grid_dims(3), 2) == 0
        boundary_idx_3 = (R3(:) == -k_grid_dims(3)/2);
        correction_factor = cos(pi * k_prime(3) * k_grid_dims(3)) ./ ...
                            exp(2i * pi * k_prime(3) * (-k_grid_dims(3)/2));
        phase_factor(boundary_idx_3) = phase_factor(boundary_idx_3) .* correction_factor;
    end
    
    exp_phase = phase_factor;
    
    % 求和: H(k') = sum_R H(R) * e^(ikR)
    H_k_prime = zeros(n_wannier, n_wannier);
    for m = 1:n_wannier
        for n = 1:n_wannier
            H_k_prime(m, n) = sum(H_fft_flat(:, m, n) .* exp_phase);
        end
    end
    
    % 对角化。WI-SCDM-f 先得到 f(E)，再通过 f^{-1} 回到真实能量。
    if use_wi_scdm_f
        H_k_prime = (H_k_prime + H_k_prime') / 2;
        evals = wi_scdm_f_inverse(eig(H_k_prime), transform_info);
    else
        evals = eig(H_k_prime);
    end
    interpolated_bands(ik, :) = sort(real(evals));
end

    efermi_eV =   efermi*27.211386246;
    interpolated_bands = interpolated_bands - efermi_eV;
    fprintf('  已平移能带，E_fermi = %.4f eV\n', efermi_eV);
fprintf('=== SCDM 插值计算完成 ===\n');
end

function options = parse_interpolation_options(varargin)
    parser = inputParser;
    addParameter(parser, 'Method', 'wi-scdm', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'TransformN', 0.5, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0);
    addParameter(parser, 'TransformA', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
        isreal(x) && isfinite(x) && x > 0));
    addParameter(parser, 'TransformEpsilon', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
        isreal(x) && isfinite(x)));
    addParameter(parser, 'TransformTopBand', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
        isreal(x) && isfinite(x) && x >= 1 && x == round(x)));
    addParameter(parser, 'InverseSamples', 20001, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && ...
        isfinite(x) && x >= 1001 && x == round(x));
    addParameter(parser, 'TransformParams', [], ...
        @(x) isempty(x) || isstruct(x));
    parse(parser, varargin{:});
    options = parser.Results;
    options.Method = normalize_interpolation_method(options.Method);
end

function method = normalize_interpolation_method(method)
    method = lower(strtrim(char(string(method))));
    switch method
        case {'wi-scdm', 'scdm', 'scdm-only', 'none'}
            method = 'wi-scdm';
        case {'wi-scdm-f', 'wiscdmf', 'wi_scdm_f', 'scdm-f'}
            method = 'wi-scdm-f';
        otherwise
            error('compute_scdm_only_bands:BadMethod', ...
                ['Unknown interpolation method "%s". Use "wi-scdm" ', ...
                 'or "wi-scdm-f".'], method);
    end
end
