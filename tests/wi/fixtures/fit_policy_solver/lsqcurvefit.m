function [x, resnorm, residual, exitflag, output] = lsqcurvefit(model, ~, E, p, varargin)
% Inject otherwise difficult-to-reach solver outcomes into writer branch tests.
global SCDM_FIT_POLICY_SOLVER
x = SCDM_FIT_POLICY_SOLVER.x;
exitflag = SCDM_FIT_POLICY_SOLVER.exitflag;
residual = model(x, E) - p;
resnorm = sum(residual.^2);
output = struct('message', 'Synthetic injected solver outcome');
end
