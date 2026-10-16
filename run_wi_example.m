function result = run_wi_example(name, stage)
%RUN_WI_EXAMPLE Run a bundled WI example without external path configuration.
% run_wi_example('analytic'); run_wi_example('si','check');
% Use one fresh MATLAB process per material stage.
if nargin < 1, name = 'analytic'; end
name = validatestring(name, {'analytic','si','soc'});
if nargin < 2, stage = 'check'; end
root = fileparts(mfilename('fullpath'));
addpath(root, '-begin'); KSSOLV_startup;
old_dir = pwd; old_path = path;
restore = onCleanup(@() restore_environment(old_dir,old_path));
if strcmp(name,'analytic')
    addpath(fullfile(root,'example','wi','analytic'),'-begin');
    result = run_analytic_demo;
else
    stage = validatestring(stage, {'check','scf','nscf','scdm','wlmf'});
    if strcmp(name,'si'), folder='si_scalar'; driver='run_si_pipeline';
    else, folder='soc_aumn'; driver='run_soc_pipeline'; end
    addpath(fullfile(root,'example','wi',folder),'-begin');
    feval(driver, stage);
    result = struct('example',name,'stage',stage,'status','completed');
end
end
function restore_environment(old_dir,old_path)
cd(old_dir); path(old_path);
end
