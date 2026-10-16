function summary = run_wi_tests
%RUN_WI_TESTS Run the bundled WI numerical and PAO-library regressions.
root = fileparts(mfilename('fullpath'));
addpath(root, '-begin'); KSSOLV_startup;
old_path = path; restore = onCleanup(@() path(old_path));
addpath(fullfile(root,'tests','wi'),'-begin');
summary = run_release_tests;
end
