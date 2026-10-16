function test_vitale_error_pipeline
% Check units, SCDM-parameter independence, provenance and cache invalidation.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'src', 'WI'), '-begin');
old_dir = pwd;
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() cd(old_dir));
cd(work_dir);

ref_ev = [-2, 0.5; -1.8, 0.8; -1.6, 1.2]; % k x band
wan_ev = ref_ev + [0.01, 0.2; -0.02, -0.3; 0.03, 0.4];
crystal = struct('efermi', 0.2, 'name', 'vitale_test', ...
    'nkpts', 1, 'nbnd', 2, 'nkxyz', [1,1,1]);
ref_ha = ref_ev.'/27.211386246 + crystal.efermi;
args = {'BandErrorMode','entangled','Method','wi-scdm'};
[mae, mx, rms, first] = calculate_err(wan_ev, ref_ha, crystal, args{:});
direct = calculate_band_error(ref_ev.', wan_ev.', 'entangled');
assert(abs(first.eta-direct) < 1e-12);
assert(abs(mae-first.mae) < 1e-12 && abs(mx-first.max_absolute_error) < 1e-12);
assert(abs(rms-first.rmse) < 1e-12);

write_scdm_input(-100, 50);
[~,~,~,second] = calculate_err(wan_ev, ref_ha, crystal, args{:});
write_scdm_input(200, 0.2);
[~,~,~,third] = calculate_err(wan_ev, ref_ha, crystal, args{:});
assert(first.eta == second.eta && second.eta == third.eta);
assert(first.eta_max == third.eta_max);
assert(has_vitale_band_error('err.txt','wi-scdm','entangled',1,0.1));
assert(~has_vitale_band_error('err.txt','wi-scdm','entangled',2,0.1));
assert(~has_vitale_band_error('err.txt','wi-scdm','entangled',1,0.2));
assert(~has_vitale_band_error('err.txt','wi-scdm-f','entangled',1,0.1));
text = fileread('err.txt');
assert(contains(text,'scdm_mu_eV 200') && contains(text,'band_error_nu_relative_fermi_eV 1'));
assert(~contains(text,'band_error_mu_eV') && ~contains(text,'band_error_sigma_eV'));

[~,~,~,custom] = calculate_err(wan_ev, ref_ev, 'Nu', 2, 'Tau', 0.2, ...
    'BandErrorMode','entangled','Method','wi-scdm-f');
assert(custom.nu == 2 && custom.tau == 0.2);
assert(has_vitale_band_error('err.txt','wi-scdm-f','entangled',2,0.2));
assert(~has_vitale_band_error('err.txt','wi-scdm-f','entangled',1,0.1));

fid = fopen('err.txt','a');
fprintf(fid,'%s\ninterpolation_method wi-scdm\nETA_entangled_eV 0.01\n',repmat('=',1,72));
fclose(fid);
assert(~has_vitale_band_error('err.txt','wi-scdm','entangled',1,0.1));
assert(has_vitale_band_error('err.txt','wi-scdm-f','entangled',2,0.2));
fprintf('PASS: Vitale unit conversion, parameter independence, provenance and cache.\n');
clear cleanup;
end

function write_scdm_input(mu,sigma)
fid = fopen('scdm.in','w');
fprintf(fid,'2 1 2\nscdm.eig erfc %.17g %.17g\n0 0 0\nKGRID 1 1 1\n',mu,sigma);
fclose(fid);
end
