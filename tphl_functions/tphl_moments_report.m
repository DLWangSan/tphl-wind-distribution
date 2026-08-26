function R = tphl_moments_report(alpha, mu, lambda, varargin)
%TPHL_MOMENTS_REPORT Raw moments, speed shape stats, power stats, optional MGF check.
%   R = TPHL_MOMENTS_REPORT(alpha, mu, lambda)
%   R = TPHL_MOMENTS_REPORT(..., 'max_r', 6, 'c_power', 1, 'v0', 0, ...
%       'mgf_check', true, 't_grid', linspace(-0.05, 0.05, 21))
%
%   Outputs struct R with fields:
%     .raw      : M1, M2, ... (vectors / named fields)
%     .speed    : mean, var, skew, kurt (Fisher-Pearson kurtosis as in paper)
%     .power    : meanP, varP, skewP, kurtP for P = c * V^3
%     .mgf      : optional table comparing E[V^k] vs numerical d^k M/dt^k|_{t=0}

    p = inputParser;
    addParameter(p, 'max_r', 6, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'c_power', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'v0', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'mgf_check', false, @islogical);
    addParameter(p, 'mgf_check_max_k', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 't_grid', linspace(-0.08, 0.08, 17), @isvector);
    parse(p, varargin{:});
    max_r = round(p.Results.max_r);
    c = p.Results.c_power;
    v0 = p.Results.v0;
    mgf_check = p.Results.mgf_check;
    mgf_check_max_k = round(p.Results.mgf_check_max_k);
    t_grid = p.Results.t_grid(:).';

    raw = struct();
    M = zeros(1, max_r);
    for r = 1:max_r
        M(r) = tphl_raw_moment(alpha, mu, lambda, r, v0);
        raw.(sprintf('M%d', r)) = M(r);
    end
    raw.M = M;

    m1 = M(1);
    m2 = M(min(2, max_r));
    speed = struct('mean', m1);
    if max_r >= 2
        speed.var = m2 - m1^2;
        speed.std = sqrt(max(speed.var, 0));
    else
        speed.var = NaN;
        speed.std = NaN;
    end
    if max_r >= 3 && speed.std > 0
        m3 = M(3);
        speed.skew = (m3 - 3*m2*m1 + 2*m1^3) / speed.std^3;
    else
        speed.skew = NaN;
    end
    if max_r >= 4 && speed.std > 0
        m4 = M(4);
        speed.kurt = (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4) / speed.std^4;
    else
        speed.kurt = NaN;
    end

    power = struct();
    need = [3, 6, 9, 12];
    for k = need
        if max_r >= k
            power.(sprintf('M%d', k)) = M(k);
        else
            power.(sprintf('M%d', k)) = tphl_raw_moment(alpha, mu, lambda, k, v0);
        end
    end
    m3p = power.M3;
    m6p = power.M6;
    m9p = power.M9;
    m12p = power.M12;
    power.meanP = c * m3p;
    power.varP = c^2 * (m6p - m3p^2);
    sdP = sqrt(max(power.varP, 0));
    if sdP > 0
        power.skewP = (c^3*m9p - 3*c^2*m6p*(c*m3p) + 2*(c*m3p)^3) / sdP^3;
        power.kurtP = (c^4*m12p - 4*c^3*m9p*(c*m3p) + 6*c^2*m6p*(c*m3p)^2 - 3*(c*m3p)^4) / sdP^4;
    else
        power.skewP = NaN;
        power.kurtP = NaN;
    end

    mgf = struct();
    if mgf_check
        kmax = min(mgf_check_max_k, max_r);
        from_mgf = nan(1, kmax);
        from_raw = M(1:kmax);
        for k = 1:kmax
            from_mgf(k) = local_mgf_derivative_at_zero(alpha, mu, lambda, k, v0);
        end
        mgf.order = (1:kmax).';
        mgf.from_raw_moment = from_raw(:);
        mgf.from_mgf_derivative = from_mgf(:);
        mgf.abs_diff = abs(from_raw(:) - from_mgf(:));
        mgf.t_grid = t_grid;
        mgf.M_t = arrayfun(@(t) tphl_mgf(alpha, mu, lambda, t, v0), t_grid);
    end

    R = struct('alpha', alpha, 'mu', mu, 'lambda', lambda, 'v0', v0, ...
        'raw', raw, 'speed', speed, 'power', power, 'mgf', mgf);
end

function dk = local_mgf_derivative_at_zero(alpha, mu, lambda, k, v0)
    f = @(t) tphl_mgf(alpha, mu, lambda, t, v0);
    h = 2e-5;
    switch k
        case 0
            dk = f(0);
        case 1
            dk = (f(h) - f(-h)) / (2*h);
        case 2
            dk = (f(h) - 2*f(0) + f(-h)) / h^2;
        case 3
            dk = (f(2*h) - 2*f(h) + 2*f(-h) - f(-2*h)) / (2*h^3);
        case 4
            dk = (f(-2*h) - 4*f(-h) + 6*f(0) - 4*f(h) + f(2*h)) / h^4;
        otherwise
            dk = local_central_deriv(f, 0, k, h);
    end
end

function d = local_central_deriv(f, t0, k, h)
    if k == 0
        d = f(t0);
        return;
    end
    d = (local_central_deriv(f, t0 + h, k - 1, h) - local_central_deriv(f, t0 - h, k - 1, h)) / (2*h);
end
