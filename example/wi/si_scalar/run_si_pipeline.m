function run_si_pipeline(step)
% One MATLAB process per step: no stale class/function state between stages.
step = validatestring(step, {'scf','nscf','scdm','wlmf','check'});
case_dir = fileparts(mfilename('fullpath'));
release_root = fileparts(fileparts(fileparts(case_dir)));
addpath(release_root, '-begin');
KSSOLV_startup;
wi_root = fullfile(release_root, 'src', 'WI');
addpath(case_dir, '-begin');
cd(case_dir);
set(groot, 'DefaultFigureVisible', 'off');
required = {'validate_vitale_error_runtime','has_vitale_band_error', ...
    'write_input_files','compute_pao_projectability','write_UNK_files', ...
    'scdm_m','compare_scdm_bands','compute_scdm_only_bands', ...
    'calculate_err','calculate_band_error','calculate_scdm_spread', ...
    'run_wannier90','write_w90_mmn','compare_w90_bands'};
for i = 1:numel(required)
    actual = which(required{i});
    expected = fullfile(wi_root, [required{i} '.m']);
    if ispc, matches = strcmpi(actual, expected); else, matches = strcmp(actual, expected); end
    assert(matches, 'Wrong function path for %s: %s (expected %s).', ...
        required{i}, actual, expected);
    fprintf('%s -> %s\n', required{i}, actual);
end
si_case_setup;
fprintf('=== %s: %s at %s ===\n', step, case_dir, ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
if strcmp(step, 'check')
    [cry, ~, metadata] = run_scf('PrepareOnly', true);
    [~, points] = si_band_path(cry.supercell);
    fprintf('Input check passed: %d atoms, %d SCF k points, %d path points.\n', ...
        size(cry.xyzlist, 1), cry.nkpts, size(points, 1));
    disp(metadata);
else
    feval(['run_' step]);
end
fprintf('=== %s finished at %s ===\n', step, ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end
