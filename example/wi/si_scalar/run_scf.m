function [cry, opt, metadata] = run_scf(varargin)
% Step 1. Diamond Si, two-atom primitive cell; no hand-tuned SCDM parameters.
parser = inputParser;
addParameter(parser, 'PrepareOnly', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
case_dir = si_case_setup;
if ~parser.Results.PrepareOnly
    diary(fullfile(case_dir, 'logs', 'run_scf.log'));
    cleanup = onCleanup(@() diary('off'));
end
kssolvpptype('default');
assert(meDef() == 1, 'This case requires KSSOLV Hartree units.');
rng(0, 'twister');
a = Atom('Si');
C = 10.263 * [0, 0.5, 0.5; 0.5, 0, 0.5; 0.5, 0.5, 0];
coefs = [0, 0, 0; 0.25, 0.25, 0.25];
cry = Crystal('supercell', C, 'atomlist', [a a], ...
    'xyzlist', coefs*C, 'ecut', 30, 'funct', 'PBE', ...
    'autokpts', [6 6 6 0 0 0], 'nbnd', 12, ...
    'smear', 'gauss', 'temperature', 300, 'name', 'Si diamond primitive');
opt = setksopt();
opt.maxscfiter = 200;
opt.scftol = 1e-8;
metadata = struct('case_name', 'si_paper_primitive', ...
    'system_name', 'Si diamond primitive', 'funct', 'PBE', ...
    'ecut', 30, 'lattice_constant_bohr', 10.263, ...
    'scf_autokpts', [6 6 6 0 0 0], ...
    'nscf_autokpts', [8 8 8 0 0 0], 'autokpts', [8 8 8 0 0 0], ...
    'nbnd_scf', 12, 'nbnd_nscf', 12, 'focc', 'erfc', ...
    'band_error_mode', 'entangled', 'nu', 1, 'tau', 0.1, ...
    'band_error_metric_definition', 'vitale2020-v1', ...
    'run_path_bands', true, 'temperature', 300, 'random_seed', 0);
metadata.target_eta_methods = {'wi-scdm', 'wi-scdm-f', 'wlmf'};
metadata.path_kpoints = si_band_path(C);
metadata.structure_source = 'Diamond primitive cell specified in this run_scf.m (a=10.263 Bohr)';
% Fail before expensive SCF if the selected PP_CHI data cannot support the run.
[~, pao_info] = compute_pao_projectability(cry, metadata.nbnd_nscf, 512, true);
metadata.pao_count = pao_info.total_count;
metadata.pao_files = pao_info.pp_files;
assert(metadata.pao_count <= metadata.nbnd_nscf, ...
    'PAO count exceeds NSCF band count; review this pseudopotential choice.');
if parser.Results.PrepareOnly, return; end
fprintf('=== Si SCF: 6x6x6, 12 bands, PBE, 30 Ha ===\n');
[cry, H, X, info] = scf(cry, opt);
converged = isfield(info, 'converge') && info.converge;
save(fullfile(case_dir, 'scf_checkpoint.mat'), 'cry', 'H', 'X', ...
    'info', 'opt', 'metadata', 'coefs', 'converged', '-v7.3');
if ~converged
    error('si:SCFNotConverged', ...
        'SCF did not converge. Check logs/run_scf.log before starting NSCF.');
end
fprintf('=== Si SCF complete ===\n');
end
