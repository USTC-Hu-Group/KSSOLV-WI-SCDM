function write_UNK_files(BX, crystal)
% write_UNK_files_new - 高效将倒空间波函数转换为UNK格式文件
% 支持自旋轨道耦合 (SOC)
%
% 逻辑：
% 1. 检查 crystal.nspin。
% 2. 如果 nspin == 1 (无自旋/共线)，生成 UNK%05d.1。
% 3. 如果 nspin == 4 (SOC)，生成 UNK%05d.NC，并处理双分量波函数。
% 
% 改动：
% 1. 适配 SOC：输入波函数行数为 2*npw。
% 2. 输出文件名适配 (.1 或 .NC)。
% 3. 实空间转换包含两个分量 (Up/Down)。
% 4. 归一化针对整个旋量 (Spinor)。

unk_dir = fullfile(pwd, 'tmp');
if ~exist(unk_dir, 'dir')
    mkdir(unk_dir);
end

% === 清理旧文件 ===
if ~isempty(dir('UNK*.1'))
    delete('UNK*.1');
    fprintf('已清理当前文件夹内的旧 UNK*.1 文件。\n');
end
if ~isempty(dir('UNK*.NC'))
    delete('UNK*.NC');
    fprintf('已清理当前文件夹内的旧 UNK*.NC 文件。\n');
end
if ~isempty(dir(fullfile(unk_dir, 'UNK*.1')))
    delete(fullfile(unk_dir, 'UNK*.1'));
    fprintf('已清理 tmp 文件夹内的旧 UNK*.1 文件。\n');
end
if ~isempty(dir(fullfile(unk_dir, 'UNK*.NC')))
    delete(fullfile(unk_dir, 'UNK*.NC'));
    fprintf('已清理 tmp 文件夹内的旧 UNK*.NC 文件。\n');
end

% 提取基本信息
if isprop(BX, 'wavefuncell') && ~isempty(BX.wavefuncell)
    nkpts = length(BX.wavefuncell);
elseif isprop(BX, 'psi')
    nkpts = 1;
else
    nkpts = 0;
end

if nkpts == 0, error('没有找到波函数数据'); end

% 检查 Spin 设置
nspin = 1;
if isprop(crystal, 'nspin')
    nspin = crystal.nspin;
end

is_soc = (nspin == 4);

if is_soc
    fprintf('检测到 SOC 模式 (nspin=4)。将生成 UNKxxxxx.NC 文件。\n');
else
    fprintf('检测到标准模式 (nspin=%d)。将生成 UNKxxxxx.1 文件。\n', nspin);
end

% 获取能带数量 (nbands)
if nkpts == 1 && isprop(BX, 'psi')
    nbands = size(BX.psi, 2);
else
    nbands = size(BX.wavefuncell{1,1}.psi, 2);
end

nx = crystal.n1;
ny = crystal.n2;
nz = crystal.n3;

fprintf('网格尺寸: %d x %d x %d\n', nx, ny, nz);
fprintf('K点数量: %d\n', nkpts);
fprintf('能带数量: %d\n', nbands);
fprintf('输出目录: %s\n', unk_dir);

% === 循环处理每个 K 点 ===
for ik = 1:nkpts
    t_start = tic;
    
    % 获取数据
    if nkpts == 1 && isprop(BX, 'psi')
        % 单 K 点 / Gamma 点模式
        psi_k = BX.psi;
        if isprop(BX, 'idxnz')
            idxnz = BX.idxnz;
        else
            error('单 K 点模式下，BX 结构体必须包含 .idxnz 字段');
        end
    else
        % 多 K 点 / wavefuncell 模式
        wave_data = BX.wavefuncell{ik, 1};
        psi_k = wave_data.psi;      % [rows, nbands]
        idxnz = wave_data.idxnz;    % [n_plane_waves, 1]
    end
    
    npw = length(idxnz);
    
    % 维度检查
    if is_soc
        if size(psi_k, 1) ~= 2 * npw
            error('k点%d (SOC模式): psi行数 (%d) 应为 idxnz长度 (%d) 的2倍', ...
                ik, size(psi_k, 1), npw);
        end
    else
        if size(psi_k, 1) ~= npw
            error('k点%d (标准模式): psi行数 (%d) 与 idxnz长度 (%d) 不一致', ...
                ik, size(psi_k, 1), npw);
        end
    end
    
    % 定义文件名
    if is_soc
        filename = sprintf('UNK%05d.NC', ik);
    else
        filename = sprintf('UNK%05d.1', ik);
    end
    
    % 打开文件
    filepath = fullfile(unk_dir, filename);
    fid = fopen(filepath, 'w');
    if fid == -1, error('无法创建文件: %s', filepath); end
    
    % 1. 写入文件头 (格式相同: nx ny nz ik nbands)
    fprintf(fid, '%5d %5d %5d %5d %5d\n', nx, ny, nz, ik, nbands);
    
    % 2. 处理所有能带
    for iband = 1:nbands
         if is_soc
            % === SOC 模式处理 (双分量) ===
        
            % --- A. 分离 Spin Up / Spin Down 分量 ---
            % 前 npw 行是 Up，后 npw 行是 Down
            psi_k_up = psi_k(1:npw, iband);
            psi_k_dw = psi_k(npw+1:2*npw, iband);
        
            % --- B. 构建倒空间网格 ---
            grid_rec_up = zeros(nx, ny, nz);
            grid_rec_dw = zeros(nx, ny, nz);
        
            grid_rec_up(idxnz) = psi_k_up;
            grid_rec_dw(idxnz) = psi_k_dw;
        
            % --- C. FFT 到实空间 ---
            psi_real_up = ifftn(grid_rec_up);
            psi_real_dw = ifftn(grid_rec_dw);
        
            % --- D. 强制实空间归一化 (针对整个 Spinor) ---
            % 先拉直为列向量，方便后续处理
            u_col = psi_real_up(:);
            d_col = psi_real_dw(:);
            
            % 计算总模长: sqrt( sum(|up|^2) + sum(|dw|^2) )
            total_norm = norm([u_col; d_col]);
        
            if total_norm > 1e-12
                u_col = u_col / total_norm;
                d_col = d_col / total_norm;
            end
        
            % --- E. 交替排列与写入 (关键修改) ---
            
            % 1. 交替混合 Up 和 Down
            % 构造矩阵 [u1, u2, ...; d1, d2, ...]
            % 展平后变成列向量: [u1, d1, u2, d2, ...]
            combined_matrix = [u_col.'; d_col.']; 
            psi_interleaved = combined_matrix(:); 
            
            % 2. 构建输出矩阵
            % 目标格式: 每一行写一个复数 (实部 虚部)
            % fprintf 是按列读取的，所以我们要构造一个 2行的矩阵：
            % Row 1: 实部
            % Row 2: 虚部
            data_to_write = [real(psi_interleaved).'; imag(psi_interleaved).'];
            
            % 3. 写入文件
            % %16.10E %16.10E\n 表示一行输出两个浮点数（实部 虚部），然后换行
            fprintf(fid, '%16.10E %16.10E\n', data_to_write);


            
        else
            % === 标准模式处理 (单分量) ===
            
            % --- A. 构建倒空间网格 ---
            psi_grid_reciprocal = zeros(nx, ny, nz);
            psi_grid_reciprocal(idxnz) = psi_k(:, iband);
            
            % --- B. FFT 到实空间 ---
            psi_grid_real = ifftn(psi_grid_reciprocal);
            
            % --- C. 强制实空间归一化 ---
            psi_grid_real = psi_grid_real / norm(psi_grid_real(:));
            
            % --- D. 向量化写入 ---
            data_to_write = [real(psi_grid_real(:))'; imag(psi_grid_real(:))'];
            
            fprintf(fid, '%16.10E %16.10E\n', data_to_write);
        end
    end
    
    fclose(fid);
    t_end = toc(t_start);
    
    % 减少打印频率
    if mod(ik, 100) == 0 || ik == 1 || ik == nkpts
        fprintf('已完成 %s (耗时 %.2fs)\n', filepath, t_end);
    end
end

fprintf('所有 UNK 文件生成完毕。\n');
end
