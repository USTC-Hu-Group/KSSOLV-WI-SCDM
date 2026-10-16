function test_scdm_fit_policy
% Characterization tests of the writer, not a material/robustness benchmark.
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_dir = pwd; old_path = path; old_warnings = warning;
cleanup = onCleanup(@() restore_environment(old_dir, old_path, old_warnings));
addpath(fullfile(repo, 'src', 'WI'), '-begin');
addpath(fullfile(repo, 'tests', 'wi', 'fixtures', 'fit_policy_pao'), '-begin');
warning('on', 'write_input_files:BadManualMu');
warning('on', 'write_input_files:BadManualSigma');
global SCDM_FIT_POLICY_FIXTURE SCDM_FIT_POLICY_SOLVER
out = fullfile(repo, 'tests', 'wi', '_work', 'fit_policy');
if ~isfolder(out), mkdir(out); end
work = tempname(out); mkdir(work); cd(work);
E = linspace(-3, 4, 8)';
p = 0.5 * erfc((E - 0.5) / 1.2);
base = struct('p', p, 'info', struct('total_count', 2, ...
    'spatial_count', 2, 'spin_factor', 1));
crystal = struct('nkpts', 1, 'nbnd', 8, 'kpts', [0 0 0], ...
    'supercell', eye(3)*8, 'xyzlist', [0 0 0], 'alist', 1, ...
    'atoms', struct('symbol', 'Si', 'anum', 14, 'venum', 4));
rows = struct('test', {}, 'expected_error', {}, 'observed_error', {}, ...
    'status', {}, 'parameter_source', {}, 'warning_identifier', {}, 'passed', {});

run_case('manual_valid', base, {2,'erfc',0,2}, '', 'manual', '');
run_case('manual_invalid_mu', base, {2,'erfc',NaN,2}, '', ...
    'manual_sanitized', 'write_input_files:BadManualMu');
run_case('manual_invalid_sigma', base, {2,'erfc',0,-1}, '', ...
    'manual_sanitized', 'write_input_files:BadManualSigma');
run_case('manual_width_below_auto_bound', base, {2,'erfc',0,0.05}, '', 'manual', '');
warning('error', 'write_input_files:BadManualMu');
run_case('strict_invalid_manual_mu', base, {2,'erfc',NaN,2}, ...
    'write_input_files:BadManualMu', 'manual', '');
warning('on', 'write_input_files:BadManualMu');
warning('error', 'write_input_files:BadManualSigma');
run_case('strict_invalid_manual_sigma', base, {2,'erfc',0,-1}, ...
    'write_input_files:BadManualSigma', 'manual', '');
warning('on', 'write_input_files:BadManualSigma');
run_case('manual_dimension_exceeds_bands', base, {9,'erfc',0,2}, ...
    'write_input_files:TooManyWannierFunctions', 'manual', '');
fixture = base; fixture.info = rmfield(fixture.info, 'total_count');
run_case('missing_pao_count', fixture, {'erfc'}, ...
    'write_input_files:MissingPAOCount', 'automatic_fit', '');
fixture = base; fixture.info.total_count = 9;
run_case('pao_count_exceeds_bands', fixture, {'erfc'}, ...
    'write_input_files:PAOCountExceedsBands', 'automatic_fit', '');
fixture = base; fixture.error_identifier = 'fit_policy:InjectedPAOFailure';
run_case('pao_provider_failure', fixture, {'erfc'}, ...
    'fit_policy:InjectedPAOFailure', 'automatic_fit', '');
fixture = base; fixture.p(4:end) = NaN;
run_case('fewer_than_four_finite_pairs', fixture, {'erfc'}, ...
    'write_input_files:TooFewFitData', 'automatic_fit', '');
fixture = base; fixture.p = linspace(0.5, 0.549, 8)';
run_case('nearly_constant_projectability', fixture, {'erfc'}, ...
    'write_input_files:FlatProjectability', 'automatic_fit', '');
fixture = base; fixture.p(1) = 1.1;
run_case('out_of_range_projectability', fixture, {'erfc'}, ...
    'write_input_files:BadProjectabilityRange', 'automatic_fit', '');
meta = run_case('automatic_fit_real_solver', base, {'erfc'}, '', 'automatic_fit', '');
assert(max(abs(meta.fit_result(:) - [0.5;1.2])) < 1e-6);
assert(meta.fit_exitflag > 0 && meta.fit_lower_bound(2) == 0.100001);
assert(max(abs(meta.fit_initial_guess(:)-[E(find(abs(p-.5)==min(abs(p-.5)),1));2])) < 1e-12);
run_case('count_only_manual_parameters', base, {'erfc',0,2}, '', 'manual', '');
assert(SCDM_FIT_POLICY_FIXTURE.observed_count_only);
run_case('isolated_count_only', base, {'isolated'}, '', 'not_applicable', '');
assert(SCDM_FIT_POLICY_FIXTURE.observed_count_only);

addpath(fullfile(repo, 'tests', 'wi', 'fixtures', 'fit_policy_solver'), '-begin');
SCDM_FIT_POLICY_SOLVER = struct('x', [0.5,1.2], 'exitflag', 0);
run_case('nonpositive_solver_exitflag', base, {'erfc'}, ...
    'write_input_files:UnreliableFit', 'automatic_fit', '');
SCDM_FIT_POLICY_SOLVER = struct('x', [0.5,0.1], 'exitflag', 1);
run_case('solver_width_at_rejection_limit', base, {'erfc'}, ...
    'write_input_files:UnreliableFit', 'automatic_fit', '');

writetable(struct2table(rows), fullfile(out, 'branch_tests.csv'));
summary = struct('test_count', numel(rows), 'passed', sum([rows.passed]), ...
    'matlab_version', version, 'work_directory', work, ...
    'scope', ['Synthetic PAO data and two injected solver outcomes; ', ...
    'one real nonlinear fit. These counts are not material success rates.']);
fid = fopen(fullfile(out, 'summary.json'), 'w');
closer = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(summary, 'PrettyPrint', true));
clear closer
fprintf('PASS: %d/%d writer policy tests. Results: %s\n', numel(rows), numel(rows), out);

    function meta = run_case(name, fixture, args, expected_error, source, expected_warning)
        folder = fullfile(work, name); mkdir(folder); cd(folder);
        % Existing outputs must not be mistaken for new results after failure.
        copyfile(fullfile(repo, 'tests', 'wi', 'test_scdm_fit_policy.m'), 'scdm.in');
        copyfile('scdm.in', 'scdm.eig');
        before = fileread('scdm.in');
        SCDM_FIT_POLICY_FIXTURE = fixture;
        lastwarn(''); observed_error = '';
        try
            write_input_files(crystal, E/27.211386246, args{:});
        catch err
            observed_error = err.identifier;
        end
        [~, warn_id] = lastwarn;
        assert(strcmp(observed_error, expected_error), '%s: unexpected exception %s', name, observed_error);
        assert(strcmp(warn_id, expected_warning), '%s: unexpected warning %s', name, warn_id);
        files = dir(fullfile(folder, 'scdm_provenance', 'run_*', 'metadata.json'));
        assert(numel(files)==1, '%s: missing or duplicate archive', name);
        meta = jsondecode(fileread(fullfile(files(1).folder, files(1).name)));
        assert(strcmp(meta.parameter_source, source));
        if isempty(expected_error)
            assert(strcmp(meta.status, 'success'));
            assert(meta.manual_parameters_sanitized == strcmp(source, 'manual_sanitized'));
            if strcmp(name, 'manual_invalid_mu'), assert(meta.mu_eV==0); end
            if strcmp(name, 'manual_invalid_sigma'), assert(meta.sigma_eV==2); end
            if strcmp(name, 'manual_width_below_auto_bound'), assert(meta.sigma_eV==0.05); end
        else
            assert(strcmp(meta.status, 'failed'));
            assert(strcmp(meta.error_identifier, expected_error));
            assert(strcmp(before, fileread('scdm.in')));
            assert(strcmp(before, fileread('scdm.eig')));
        end
        rows(end+1) = struct('test', name, 'expected_error', expected_error, ...
            'observed_error', observed_error, 'status', meta.status, ...
            'parameter_source', source, 'warning_identifier', warn_id, 'passed', true);
        cd(work);
    end
end

function restore_environment(old_dir, old_path, old_warnings)
cd(old_dir); path(old_path); warning(old_warnings);
clear global SCDM_FIT_POLICY_FIXTURE SCDM_FIT_POLICY_SOLVER
end
