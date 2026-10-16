function [energies_ev, diagnostics] = wi_scdm_f_inverse(transformed_ev, params)
%WI_SCDM_F_INVERSE Recover physical energies from WI-SCDM-f eigenvalues.
%
% Python-compatible default: invert transform.py's fun(a,n,s,x) with the
% step-limited Newton iteration, checking residuals before updating. Values above zero
% are clipped to epsilon because f is flat for E >= epsilon.
% The optional diagnostics describe residuals and clipping, not physical accuracy.

    validate_inverse_params(params);
    transformed_ev = double(real(transformed_ev));

    a = params.a;
    epsilon = params.epsilon;

    target = transformed_ev;
    positive_mask = target > 0;
    if any(positive_mask(:))
        max_positive = max(target(positive_mask));
        if max_positive > 1e-8 * max(a, 1)
            warning('wi_scdm_f_inverse:PositiveTransformedEnergy', ...
                ['Encountered transformed eigenvalues above zero ', ...
                 '(max %.6g eV). Clipping them to epsilon.'], ...
                max_positive);
        end
        target(positive_mask) = 0;
    end

    dxmax = a / 10;
    resmax = 1e-14;
    energies_ev = target + epsilon;
    residual = wi_scdm_f_transform(energies_ev, params) - target;
    niter = 0;

    while max(abs(residual(:))) > resmax && niter < 50
        derivative = wi_scdm_f_derivative(energies_ev, params);
        step = zeros(size(target));
        mask = abs(residual) > resmax & abs(derivative) ~= 0;
        step(mask) = -residual(mask) ./ derivative(mask);
        step(step > dxmax) = dxmax;
        step(step < -dxmax) = -dxmax;
        energies_ev = energies_ev + step;
        niter = niter + 1;
        residual = wi_scdm_f_transform(energies_ev, params) - target;
    end

    max_residual = max(abs(residual(:)));
    diagnostics = struct('iterations', niter, 'residual_tolerance_eV', resmax, ...
        'max_residual_eV', max_residual, 'converged', max_residual <= resmax, ...
        'clipped_positive_count', nnz(positive_mask), ...
        'zero_target_count', nnz(target == 0), ...
        'clipping_warning_threshold_eV', 1e-8 * max(a,1));
    if max_residual > resmax
        warning('wi_scdm_f_inverse:NewtonNotConverged', ...
            ['Newton inverse did not converge to %.1e after 50 ', ...
             'iterations; max residual = %.6e eV.'], ...
            resmax, max_residual);
    end
end

function derivative = wi_scdm_f_derivative(energies_ev, params)
    a = params.a;
    n = params.n;
    epsilon = params.epsilon;
    erf_half = erf(n / 2);
    if abs(erf_half) < eps
        error('wi_scdm_f_inverse:BadN', ...
            'TransformN is too small; erf(n/2) is numerically zero.');
    end

    y = energies_ev - epsilon;
    derivative = zeros(size(y));
    low_mask = y <= -a;
    derivative(low_mask) = 1;

    mid_mask = y > -a & y < 0;
    if any(mid_mask(:))
        ym = y(mid_mask);
        derivative(mid_mask) = 0.5 - ...
            erf(n * (0.5 + ym / a)) ./ (2 * erf_half);
    end
end

function validate_inverse_params(params)
    required = {'a', 'n', 'epsilon'};
    for index = 1:numel(required)
        if ~isstruct(params) || ~isfield(params, required{index})
            error('wi_scdm_f_inverse:MissingParam', ...
                'Transform params must contain field %s.', required{index});
        end
    end
    if isfield(params, 'inverse_samples') && ~isempty(params.inverse_samples) && ...
            (~isnumeric(params.inverse_samples) || ...
            ~isscalar(params.inverse_samples) || ...
            ~isfinite(params.inverse_samples) || params.inverse_samples < 1001)
        error('wi_scdm_f_inverse:BadInverseSamples', ...
            'inverse_samples must be at least 1001.');
    end
end
