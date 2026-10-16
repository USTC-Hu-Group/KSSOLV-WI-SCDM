function [p, info] = compute_pao_projectability(~, ~, ~, count_only)
% Synthetic test double, never a source of manuscript material results.
global SCDM_FIT_POLICY_FIXTURE
if isfield(SCDM_FIT_POLICY_FIXTURE, 'error_identifier')
    error(SCDM_FIT_POLICY_FIXTURE.error_identifier, 'Injected PAO failure.');
end
SCDM_FIT_POLICY_FIXTURE.observed_count_only = count_only;
p = SCDM_FIT_POLICY_FIXTURE.p;
info = SCDM_FIT_POLICY_FIXTURE.info;
if count_only, p = []; end
end
