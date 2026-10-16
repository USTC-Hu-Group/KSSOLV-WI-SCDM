function case_dir = si_case_setup
% Anchor every checkpoint to this case, never to another MATLAB path entry.
case_dir = fileparts(mfilename('fullpath'));
cd(case_dir);
for name = {'logs', 'results', 'tmp', 'wannier90_post'}
    folder = fullfile(case_dir, name{1});
    if ~isfolder(folder), mkdir(folder); end
end
validate_vitale_error_runtime;
end
