function run_soc_pipeline(step)
%RUN_SOC_PIPELINE Execute one existing AuMn stage in a fresh MATLAB process.
step = validatestring(step, {'check','scf','nscf','scdm','wlmf'});
case_dir = fileparts(mfilename('fullpath'));
release_root = fileparts(fileparts(fileparts(case_dir)));
addpath(release_root, '-begin');
KSSOLV_startup;
addpath(case_dir, '-begin');
cd(case_dir);
set(groot, 'DefaultFigureVisible', 'off');
% Must precede all BlochHam construction and checkpoint loading.
enable_soc_blochham_compat;
validate_vitale_error_runtime;
if strcmp(step, 'check')
    run_scf('PrepareOnly', true);
else
    feval(['run_' step]);
end
fprintf('SOC stage %s completed.\n', step);
end
