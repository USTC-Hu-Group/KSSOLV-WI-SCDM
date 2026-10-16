function test_vitale_error_summary
% Exercise actual evaluator records through the population/pair collector.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'src', 'WI'), '-begin');
old_dir = pwd;
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() cd(old_dir));
case_dir = fullfile(work_dir,'case_a');
mkdir(case_dir);
cd(case_dir);
ref = [-2; -1; 0];
[~,~,~,scdm] = calculate_err(ref+0.2, ref, 'BandErrorMode','entangled', ...
    'Method','wi-scdm');
[~,~,~,scdm_f] = calculate_err(ref+0.3, ref, 'BandErrorMode','entangled', ...
    'Method','wi-scdm-f');
[~,~,~,mlwf] = calculate_err(ref+0.1, ref, 'BandErrorMode','entangled', ...
    'Method','wlmf');
% A more recent alternate-window record must not replace the default result.
calculate_err(ref+2, ref, 'BandErrorMode','entangled', 'Nu',2,'Method','wlmf');
fid = fopen('err.txt','a');
fprintf(fid,'%s\ninterpolation_method wlmf\nETA_entangled_eV 0.9\n',repmat('=',1,72));
fclose(fid);
[summary, paired] = collect_eta_summary(work_dir);
assert(height(summary) == 3 && height(paired) == 1);
assert(abs(paired.wlmf_eta_eV-mlwf.eta) < 1e-11);
assert(abs(paired.delta_eta_wlmf_minus_wi_scdm_eV-(mlwf.eta-scdm.eta)) < 1e-11);
assert(abs(paired.delta_eta_wlmf_minus_wi_scdm_f_eV-(mlwf.eta-scdm_f.eta)) < 1e-11);
assert(paired.delta_eta_wlmf_minus_wi_scdm_eV < 0);
assert(isfile(fullfile(work_dir,'eta_summary.csv')));
assert(isfile(fullfile(work_dir,'eta_pair_summary.csv')));
empty_dir = fullfile(work_dir,'empty');
mkdir(empty_dir);
[empty_summary, empty_pairs] = collect_eta_summary(empty_dir);
assert(height(empty_summary) == 0 && height(empty_pairs) == 0);
fprintf('PASS: Vitale summary excludes old/other-window metrics and pairs deltas.\n');
clear cleanup;
end
