function [cry, opt, metadata] = run_scf(varargin)
% Atom(symbol, moment, theta_deg, phi_deg); angles are input in degrees.
parser = inputParser;
addParameter(parser, 'PrepareOnly', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
if ~isfolder('logs'), mkdir('logs'); end
diary(fullfile('logs', 'run_scf.log'));
cleanup = onCleanup(@() diary('off'));
validate_vitale_error_runtime;
enable_soc_blochham_compat;
kssolvpptype('default'); % lspinorb selects ppdata/default_FR in PpData.
assert(meDef() == 1, 'This input assumes the Hartree KSSOLV convention.');

a_Au = Atom('Au', 0, 0, 0);
a_Mn = Atom('Mn', 3, 0, 0);
pp_Au = PpData(a_Au, true);
assert(pp_Au.info.has_so && strcmpi(pp_Au.info.relativistic, 'full'), ...
    'SOC requires fully relativistic Au pseudopotentials.');
assert(abs(pp_Au.venum - 19) < 1e-10, ...
    'Au pseudopotential valence differs from the generation manifest.');
pp_Mn = PpData(a_Mn, true);
assert(pp_Mn.info.has_so && strcmpi(pp_Mn.info.relativistic, 'full'), ...
    'SOC requires fully relativistic Mn pseudopotentials.');
assert(abs(pp_Mn.venum - 15) < 1e-10, ...
    'Mn pseudopotential valence differs from the generation manifest.');
atomlist = [a_Mn a_Au];
C = [
          5.795469880811       0.000000000000       0.000000000000
          0.000000000000       5.795469880811       0.000000000000
          0.000000000000       0.000000000000       5.795469880811
];
coefs = [
          0.500000000000       0.500000000000       0.500000000000
          0.000000000000       0.000000000000       0.000000000000
];
xyzlist = coefs * C; % XSF lattice/coordinates converted from Angstrom to Bohr.
scf_autokpts = [4 4 4 0 0 0];
nscf_autokpts = [6 6 6 0 0 0];
nbnd_scf = 42;
nbnd_nscf = 90;
temperature = 300; % K; finite Gaussian smearing also supports metallic cases.
cry = Crystal('nspin', 4, 'noncolin', true, 'lspinorb', true, ...
    'domag', 1, 'funct', 'PBE', ...
    'supercell', C, 'atomlist', atomlist, 'xyzlist', xyzlist, ...
    'nel', 34, 'nbnd', nbnd_scf, ...
    'ecut', 40, 'name', 'AuMn', ...
    'autokpts', scf_autokpts, 'smear', 'gauss', 'temperature', temperature);
validate_soc_kssolv_runtime(cry);

opt = setksopt();
opt.maxscfiter = 200;
opt.betamix = 0.3;
opt.mixdim = 8;
opt.scftol = 1e-8;

metadata = struct();
metadata.case_name = 'soc_aumn';
metadata.system_name = 'AuMn';
metadata.computational_setting = 'SOC-included calculations';
metadata.source_kind = 'xsf';
metadata.source_xsf = 'structure.xsf';
metadata.pptype = 'default';
metadata.pseudo_subdir = 'default_FR';
metadata.pseudo_files = struct();
metadata.pseudo_files.Au = 'Au_ONCV_PBE_FR-1.3.upf';
metadata.pseudo_files.Mn = 'Mn_ONCV_PBE_FR-1.3.upf';
metadata.nspin = 4;
metadata.noncolin = true;
metadata.lspinorb = true;
metadata.domag = 1;
metadata.ecut = 40;
metadata.temperature = temperature;
metadata.scf_autokpts = scf_autokpts;
metadata.nscf_autokpts = nscf_autokpts;
metadata.autokpts = nscf_autokpts;
metadata.nbnd_scf = nbnd_scf;
metadata.nbnd_nscf = nbnd_nscf;
metadata.estimated_num_wann_from_pao = 72;
metadata.focc = 'erfc';
metadata.band_error_mode = 'entangled';
metadata.band_error_metric_definition = 'vitale2020-v1';
metadata.target_eta_methods = {'wi-scdm', 'wi-scdm-f', 'wlmf'};
% A common comparison path in THIS XSF reciprocal cell; K1..K4 are generic
% labels, not an automatic crystallographic high-symmetry assignment.
metadata.path_kpoints = {
    0.0, 0.0, 0.0, 25, 'G';
    0.5, 0.0, 0.0, 25, 'K1';
    0.5, 0.5, 0.0, 25, 'K2';
    0.5, 0.5, 0.5, 25, 'K3';
    0.0, 0.0, 0.0, 25, 'G';
    0.0, 0.0, 0.5, 25, 'K4';
};
metadata.run_path_bands = true;
metadata.plot_cubes = false;
metadata.plot_locality = false;

% Count only: no UNK files or expensive projectability fit before SCF.
[~, pao_info] = compute_pao_projectability(cry, nbnd_nscf, cry.nkpts, true);
assert(pao_info.spin_factor == 2, 'Expected two PAO spin components.');
assert(pao_info.total_count == metadata.estimated_num_wann_from_pao, ...
    'Runtime PAO count differs from the generation manifest.');
assert(nbnd_nscf >= pao_info.total_count, 'NSCF nbnd must cover automatic num_wann.');
metadata.pao_preflight = pao_info;
if parser.Results.PrepareOnly
    fprintf('SOC input checks passed for %s; SCF not executed.\n', metadata.case_name);
    return;
end

[cry, H, X, info] = scf(cry, opt);
save('scf_checkpoint.mat', 'cry', 'H', 'X', 'info', 'opt', ...
    'metadata', 'coefs', '-v7.3');
if ~isfield(info, 'converge') || ~info.converge
    error('run_scf:NotConverged', ...
        'SCF did not converge. Checkpoint saved; inspect before running NSCF.');
end
fprintf('SCF checkpoint saved for %s.\n', metadata.case_name);
end
