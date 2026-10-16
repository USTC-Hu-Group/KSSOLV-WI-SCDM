function test_calculate_err_trims_reference_bands
%TEST_CALCULATE_ERR_TRIMS_REFERENCE_BANDS Compare SCDM bands with a larger DFT set.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'src', 'WI'), '-begin');

nk = 21;
nwan = 8;
nbnd = 12;
hartree_to_ev = 27.211386246;

bands_scdm = reshape(1:(nk * nwan), nk, nwan) / 100;
extra_dft_bands = 10 + reshape(1:(nk * (nbnd - nwan)), nk, nbnd - nwan) / 100;
dft_ev = [bands_scdm, extra_dft_bands];
energyBands = dft_ev.' / hartree_to_ev;

crystal = struct('efermi', 0, 'name', 'unit_test', ...
    'nkpts', nk, 'nbnd', nbnd, 'nkxyz', [nk 1 1]);

old_dir = pwd;
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() cd(old_dir));
cd(work_dir);

[err_aver, err_max, err_rms] = calculate_err( ...
    bands_scdm, energyBands, crystal, 'Method', 'wi-scdm-f');

assert(err_aver < 1e-12);
assert(err_max < 1e-12);
assert(err_rms < 1e-12);
assert(isfile('err.txt'));
err_text = fileread('err.txt');
assert(contains(err_text, 'interpolation_method wi-scdm-f'));
assert(contains(err_text, 'reference_dft_band_indices 1 2 3 4 5 6 7 8'));

fprintf('PASS: calculate_err trims extra DFT bands to the SCDM band count.\n');
clear cleanup
end
