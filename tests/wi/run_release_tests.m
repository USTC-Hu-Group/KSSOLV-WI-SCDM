function summary = run_release_tests
%RUN_RELEASE_TESTS WI regressions, including PAO lookup in the bundled KSSOLV.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path = path; old_dir = pwd; old_visible = get(groot,'DefaultFigureVisible');
restore = onCleanup(@() restore_environment(old_path,old_dir,old_visible));
addpath(root, '-begin'); KSSOLV_startup;
addpath(fullfile(root,'tests','wi'), '-begin');
set(groot,'DefaultFigureVisible','off');
work=fullfile(root,'tests','wi','_work'); if ~isfolder(work), mkdir(work); end
names = {'test_scdm_shifted_kgrid','test_calculate_band_error_modes', ...
    'test_calculate_err_trims_reference_bands','test_vitale_error_pipeline', ...
    'test_vitale_error_summary','test_scdm_provenance','audit_wi_scdm_f_inverse'};
passed = {}; skipped = {};
for i=1:numel(names)
    fprintf('\nRunning %s\n',names{i});
    feval(names{i}); passed{end+1}=names{i}; %#ok<AGROW>
    close all;
end
results = runtests(fullfile(root,'tests','wi','test_write_w90_mmn_from_unk.m'));
assert(all([results.Passed]),'release:MMNTestFailure','One or more overlap tests failed.');
passed{end+1}='test_write_w90_mmn_from_unk (3 tests)';
if license('test','Optimization_Toolbox') && exist('lsqcurvefit','file')==2
    test_scdm_fit_policy;
    passed{end+1}='test_scdm_fit_policy';
else
    skipped{end+1}='test_scdm_fit_policy: Optimization Toolbox unavailable';
end
test_compute_pao_projectability_ppvar_funct;
passed{end+1}='test_compute_pao_projectability_ppvar_funct';
summary=struct('matlab_version',version,'passed',{passed},'skipped',{skipped});
save(fullfile(work,'test_summary.mat'),'summary');
fprintf('\nRelease regressions passed: %d groups; skipped: %d.\n',numel(passed),numel(skipped));
end
function restore_environment(old_path,old_dir,old_visible)
cd(old_dir); path(old_path); set(groot,'DefaultFigureVisible',old_visible);
end
