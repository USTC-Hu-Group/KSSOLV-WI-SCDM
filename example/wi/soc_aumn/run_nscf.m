function [crystal, options_nscf, metadata] = run_nscf(varargin)
parser = inputParser;
addParameter(parser, 'PrepareOnly', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
enable_soc_blochham_compat;
kssolvpptype('default');
ensure_dir('logs');
ensure_dir('tmp');
diary(fullfile('logs', 'run_nscf.log'));
cleanup = onCleanup(@() diary('off'));

load('scf_checkpoint.mat', 'cry', 'H', 'info', 'metadata');
fprintf('=== NSCF start: %s ===\n', metadata.case_name);

assert(isfield(info, 'converge') && info.converge, 'SCF checkpoint is not converged.');
assert(cry.nspin == 4 && cry.noncolin && cry.lspinorb, 'Expected an SOC SCF checkpoint.');
validate_soc_kssolv_runtime(cry);
% autokpts creates FRACTIONAL k points. Convert exactly once with finalize.
crystal = set(cry, 'scfkpts', cry.kpts, ...
    'autokpts', metadata.nscf_autokpts, 'nbnd', metadata.nbnd_nscf);
crystal = finalize(crystal);
crystal.ppvar = []; % Rebuild radial tables for the NEW k grid.
assert(crystal.nkpts == prod(metadata.nscf_autokpts(1:3)));
assert(crystal.nbnd >= metadata.estimated_num_wann_from_pao);

rho = H.rho;
scf_info = info;
options_nscf = setksopt();
options_nscf.rho0 = rho;
options_nscf.X0 = [];
options_nscf.scfX = [];
options_nscf.iterscf = length(scf_info.Etotvec);

assert(iscell(rho) && numel(rho) == 4, 'SOC NSCF requires all four SCF density components.');
if parser.Results.PrepareOnly
    fprintf('SOC NSCF setup checked; eigensolver not executed.\n');
    return;
end
[crystal, BH, BX, info] = nscf4c(crystal, options_nscf);
write_UNK_files(BX, crystal);

save('nscf_checkpoint.mat', 'crystal', 'BH', 'BX', 'info', ...
    'options_nscf', 'rho', 'metadata', '-v7.3');
fprintf('NSCF checkpoint saved: %s\n', fullfile(pwd, 'nscf_checkpoint.mat'));

if isfield(metadata, 'run_path_bands') && metadata.run_path_bands
    try
        [crystal_path, energyBands] = run_path_bands_with_refreshed_ppvar( ...
            crystal, options_nscf, metadata.path_kpoints, rho);
        save('nscf_path_checkpoint.mat', 'crystal_path', ...
            'energyBands', 'metadata', '-v7.3');
        fprintf('Path-band checkpoint saved: %s\n', ...
            fullfile(pwd, 'nscf_path_checkpoint.mat'));
    catch path_error
        warning('run_nscf:PathBandsFailed', ...
            'Path-band NSCF failed: %s', path_error.message);
    end
end

fprintf('=== NSCF done: %s ===\n', metadata.case_name);
end

function [crystal_path, energyBands] = run_path_bands_with_refreshed_ppvar( ...
        crystal, options_nscf, path_kpoints, rho)
symkpts = reshape([path_kpoints{:, 1:3}], [], 3);
nkpts = [path_kpoints{:, 4}];
kpoints = [path_kpoints{1, 1:3}];
for i = 1:size(path_kpoints, 1) - 1
    nk = nkpts(i) - 1;
    dk = (symkpts(i + 1, 1:3) - symkpts(i, 1:3)) .* ((1:nk)' / nk);
    kpoints = [kpoints; symkpts(i, 1:3) + dk]; %#ok<AGROW>
end

options_nscf.rho0 = rho;
options_nscf.scfX = [];
options_nscf.X0 = [];
crystal_path = set(crystal, 'kpts', kpoints, 'scfkpts', crystal.kpts);
crystal_path = finalize(crystal_path);

% eband() reuses crystal.ppvar after changing k points.  For some ionic
% crystals this leaves the nonlocal pseudopotential interpolation table too
% short near the path boundary and getvnl can index past the table.  Clearing
% ppvar forces iterinit/nscf4c to rebuild it for the actual path k points.
crystal_path.ppvar = [];

[crystal_path, ~, ~, info_path] = nscf4c(crystal_path, options_nscf);
energyBands = reshape(info_path.Eigvals, crystal_path.nbnd, []);
end

function ensure_dir(dirname)
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
end
