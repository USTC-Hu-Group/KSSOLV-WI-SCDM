function transformed_ev = wi_scdm_f_transform(energies_ev, params)
%WI_SCDM_F_TRANSFORM Apply the Wu et al. WI-SCDM-f energy transform.
%
% The transform is defined so that f(epsilon) = 0. Energies above epsilon
% are flattened to zero, energies below epsilon-a are shifted linearly, and
% the interval [epsilon-a, epsilon] is smoothed by an erf transition.

    validate_transform_params(params);
    energies_ev = double(energies_ev);

    a = params.a;
    n = params.n;
    epsilon = params.epsilon;
    erf_half = erf(n / 2);
    if abs(erf_half) < eps
        error('wi_scdm_f_transform:BadN', ...
            'TransformN is too small; erf(n/2) is numerically zero.');
    end

    y = energies_ev - epsilon;
    transformed_ev = zeros(size(y));

    low_mask = y < -a;
    transformed_ev(low_mask) = y(low_mask) + a / 2;

    mid_mask = ~low_mask & y < 0;
    if any(mid_mask(:))
        ym = y(mid_mask);
        u = n * (0.5 + ym / a);
        u0 = n / 2;
        transformed_ev(mid_mask) = ym / 2 ...
            - a ./ (2 * erf_half * n) .* ( ...
            u .* erf(u) - u0 * erf_half + ...
            (exp(-u.^2) - exp(-u0^2)) / sqrt(pi));
    end
end

function validate_transform_params(params)
    required = {'a', 'n', 'epsilon'};
    for index = 1:numel(required)
        if ~isstruct(params) || ~isfield(params, required{index})
            error('wi_scdm_f_transform:MissingParam', ...
                'Transform params must contain field %s.', required{index});
        end
    end
    if ~isnumeric(params.a) || ~isscalar(params.a) || ...
            ~isreal(params.a) || ~isfinite(params.a) || params.a <= 0
        error('wi_scdm_f_transform:BadA', ...
            'Transform parameter a must be a positive finite scalar.');
    end
    if ~isnumeric(params.n) || ~isscalar(params.n) || ...
            ~isreal(params.n) || ~isfinite(params.n) || params.n <= 0
        error('wi_scdm_f_transform:BadN', ...
            'Transform parameter n must be a positive finite scalar.');
    end
    if ~isnumeric(params.epsilon) || ~isscalar(params.epsilon) || ...
            ~isreal(params.epsilon) || ~isfinite(params.epsilon)
        error('wi_scdm_f_transform:BadEpsilon', ...
            'Transform parameter epsilon must be a finite scalar.');
    end
end
