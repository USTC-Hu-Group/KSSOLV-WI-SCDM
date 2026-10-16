function test_calculate_band_error_modes
% Numerical regression against Eqs. (18), (19), (24)-(26), not erfc weights.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'src', 'WI'), '-begin');

% Printed Eq. (18) is a norm; ordinary RMS is a separate statistic.
[eta, d] = calculate_band_error([-10, -10], [-7, -6], 'isolated');
assert(abs(eta - 5) < 1e-12);
assert(abs(d.rmse - sqrt(12.5)) < 1e-12);
assert(d.eta_max == 4 && d.mae == 3.5 && d.max_absolute_error == 4);
[eta2, d2] = calculate_band_error([-10, -10, -10, -10], ...
    [-7, -6, -7, -6], 'isolated');
assert(abs(eta2 - sqrt(2)*eta) < 1e-12 && abs(d2.rmse-d.rmse) < 1e-12);

% Known FD occupations: [1/2,3/4] vs [1/4,1/4]. Unequal errors distinguish
% geometric-mean weights from one-sided weights, and w*delta from sqrt(w)*delta.
a = log(3);
ref = [0; -a];
wan = [a; a];
w = [sqrt(1/8); sqrt(3)/4];
[eta, d] = calculate_band_error(ref, wan, 'entangled', 'Nu', 0, 'Tau', 1);
expected = a * sqrt((w(1) + 4*w(2)) / sum(w));
assert(abs(eta - expected) < 1e-12);
assert(max(abs(d.weight_matrix - w)) < 1e-12);
assert(abs(d.eta_max - 2*a*w(2)) < 1e-12);
assert(max(abs(d.band_weighted_error - [a; 2*a])) < 1e-12);
assert(d.max_absolute_error > d.eta_max);

[swapped, ds] = calculate_band_error(wan, ref, 'entangled', 'Nu', 0, 'Tau', 1);
[shifted, dt] = calculate_band_error(ref+12, wan+12, 'entangled', 'Nu', 12, 'Tau', 1);
assert(abs(swapped-eta) < 1e-12 && abs(ds.eta_max-d.eta_max) < 1e-12);
assert(abs(shifted-eta) < 1e-12 && max(abs(dt.weight_matrix-w)) < 1e-12);

[default_eta, dd] = calculate_band_error(ref, wan, 'entangled');
explicit_eta = calculate_band_error(ref, wan, 'entangled', 'Nu', 1, 'Tau', 0.1);
assert(default_eta == explicit_eta && dd.nu == 1 && dd.tau == 0.1);
assert(strcmp(dd.metric_definition, 'vitale2020-v1'));
assert(~isfield(dd, 'mu') && ~isfield(dd, 'sigma'));

% The FD product underflows here, but its geometric mean must remain positive.
[tail_eta, tail] = calculate_band_error([0; 40], [0.01; 41], 'entangled');
assert(isfinite(tail_eta) && tail.weight_matrix(2) > 0);
assert(abs(tail.band_weighted_error(2)-1) < 1e-12);
assert_throws(@() calculate_band_error(1e6, 1e6+1, 'entangled'), ...
    'calculate_band_error:ZeroWeight');
assert_throws(@() calculate_band_error(ref, wan, 'entangled', 'Mu', 1, 'Sigma', 2), ...
    'calculate_band_error:SCDMParametersNotEvaluationParameters');
assert_throws(@() calculate_band_error(ref, wan, 'entangled', 'Tau', 0), '');
assert_throws(@() calculate_band_error(ref, wan, 'entangled', 'Nu', NaN), '');

plot_prefix = fullfile(tempname, 'vitale_metric');
[~, dp] = calculate_band_error(ref, wan, 'entangled', 'plot_band_error', true, ...
    'OutputPrefix', plot_prefix, 'Visible', 'off');
assert(isfile(dp.plot_files.band_error_png) && isfile(dp.plot_files.energy_weight_png));
fprintf('PASS: Vitale formulas, symmetry, defaults, tails, validation and plots.\n');
end

function assert_throws(callback, identifier)
thrown = false;
try
    callback();
catch err
    thrown = true;
    if ~isempty(identifier)
        assert(strcmp(err.identifier, identifier), 'Unexpected error: %s', err.message);
    end
end
assert(thrown, 'Expected error was not raised.');
end
