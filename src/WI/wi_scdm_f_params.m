function params = wi_scdm_f_params(eigenvalues_ev, varargin)
%WI_SCDM_F_PARAMS Choose parameters for the WI-SCDM-f transform.
%
% params = wi_scdm_f_params(eigenvalues_ev)
% params = wi_scdm_f_params(eigenvalues_ev, 'TransformN', 0.5, ...)
%
% eigenvalues_ev is [nband x nk] in eV. In the Python-compatible default,
% n = 0.5. If no explicit epsilon is supplied, epsilon is the maximum
% energy of the chosen top band, matching transform.py's fallback
% s = max(prefix.eig). compare_scdm_bands supplies the path-band maximum
% when it is available, matching transform.py's s = max(aiida_band.dat).
% a = 4 * (max_k E_top(k) - min_k E_top(k)) unless manually supplied.

    validateattributes(eigenvalues_ev, {'numeric'}, ...
        {'2d', 'nonempty', 'real', 'finite'}, ...
        'wi_scdm_f_params', 'eigenvalues_ev');

    [nband, ~] = size(eigenvalues_ev);

    parser = inputParser;
    addParameter(parser, 'TransformN', 0.5, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0);
    addParameter(parser, 'TransformA', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
        isreal(x) && isfinite(x) && x > 0));
    addParameter(parser, 'TransformEpsilon', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
        isreal(x) && isfinite(x)));
    addParameter(parser, 'TransformTopBand', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
        isreal(x) && isfinite(x) && x >= 1 && x == round(x)));
    addParameter(parser, 'InverseSamples', 20001, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && ...
        isfinite(x) && x >= 1001 && x == round(x));
    parse(parser, varargin{:});
    options = parser.Results;

    top_band = options.TransformTopBand;
    if isempty(top_band)
        top_band = nband;
    end
    top_band = double(top_band);
    if top_band > nband
        error('wi_scdm_f_params:BadTopBand', ...
            'TransformTopBand=%d exceeds nband=%d.', top_band, nband);
    end

    top_energies = eigenvalues_ev(top_band, :);

    epsilon = options.TransformEpsilon;
    if isempty(epsilon)
        epsilon = max(top_energies);
    end

    a = options.TransformA;
    if isempty(a)
        top_span = max(top_energies) - min(top_energies);
        a = 4 * top_span;
        if ~(isfinite(a) && a > 0)
            full_span = max(eigenvalues_ev(:)) - min(eigenvalues_ev(:));
            a = max(0.1, 0.05 * max(full_span, 1.0));
            warning('wi_scdm_f_params:FlatTopBand', ...
                ['Top-band energy span is zero or invalid; using ', ...
                 'fallback TransformA = %.6g eV.'], a);
        end
    end

    params = struct();
    params.method = 'wi-scdm-f';
    params.n = double(options.TransformN);
    params.a = double(a);
    params.epsilon = double(epsilon);
    params.top_band = top_band;
    params.inverse_samples = double(options.InverseSamples);
    params.inverse_method = 'newton';
    params.compatibility = 'transform.py';
    params.energy_unit = 'eV';
end
