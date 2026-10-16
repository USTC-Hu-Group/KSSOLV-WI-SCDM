function test_compute_pao_projectability_ppvar_funct
%TEST_COMPUTE_PAO_PROJECTABILITY_PPVAR_FUNCT Verify PAO lookup follows the loaded PP functional.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'src', 'WI'), '-begin');

addpath(repo_root, '-begin');
KSSOLV_startup;
kssolvpptype('default');

crystal = struct();
crystal.atoms = struct('symbol', 'Si');
crystal.alist = 1;
crystal.xyzlist = [0 0 0];
crystal.funct = 'PZ';
crystal.ppvar = struct('funct', 'PBE');
crystal.nspin = 1;

[proj_vals, pao_info] = compute_pao_projectability(crystal, 12, 1, true);

assert(isempty(proj_vals));
assert(isfield(pao_info, 'total_count') && pao_info.total_count > 0);
assert(contains(pao_info.pp_library, [filesep 'default']));
assert(contains(lower(pao_info.pp_files{1}), 'si_oncv_pbe'));

fprintf('PASS: PZ crystal.funct can use the actually loaded PBE default PAO.\n');
end
