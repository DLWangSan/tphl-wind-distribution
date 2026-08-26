function Mr = tphl_raw_moment(alpha, mu, lambda, r, v0)
%TPHL_RAW_MOMENT Raw moment E[V^r] under truncated TPHL (V >= v0).
%   Mr = TPHL_RAW_MOMENT(alpha, mu, lambda, r, v0)
%
%   Implements Supplementary Eq. (S8):
%     E[V^r | V >= v0] = (1/S0) * int_{v0}^infty v^r f_PHL(v) dv
%   equivalently int v^r f_TPHL(v) dv.
%
%   Reference: Wang et al. (2026) OLAR 0163, Main Eq. (7), Supp. (S8).

    if nargin < 5, v0 = 0; end
    if r <= 0 || ~isfinite(r)
        error('r must be a positive finite number.');
    end
    if alpha <= 0 || lambda <= 0
        error('alpha and lambda must be positive.');
    end

    v99 = tphl_quantile(0.999, alpha, mu, lambda, v0);
    vmax = max(v99, v0) + 10 / max(lambda, 1e-3);
    fv = @(v) tphl_pdf(v, alpha, mu, lambda, v0);

    Mr = integral(@(v) (v.^r) .* fv(v), v0, vmax, ...
        'RelTol', 1e-8, 'AbsTol', 1e-11);
end
