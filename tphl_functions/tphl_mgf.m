function M = tphl_mgf(alpha, mu, lambda, t, v0)
%TPHL_MGF Moment generating function M_V(t) = E[exp(t V)] under TPHL.
%   M = TPHL_MGF(alpha, mu, lambda, t, v0)
%
%   Implements Supplementary Eq. (S9):
%     M_V(t) = (1/S0) * int_{v0}^infty exp(t v) f_PHL(v) dv
%
%   For |t| large, the integral may diverge; caller should keep t in a
%   neighborhood of 0 where M_V(t) exists.
%
%   Reference: Wang et al. (2026) OLAR 0163, Supp. (S9).

    if nargin < 5, v0 = 0; end
    if alpha <= 0 || lambda <= 0
        error('alpha and lambda must be positive.');
    end

    v99 = tphl_quantile(0.999, alpha, mu, lambda, v0);
    vmax = max(v99, v0) + 10 / max(lambda, 1e-3);
    fv = @(v) tphl_pdf(v, alpha, mu, lambda, v0);

    M = integral(@(v) exp(t .* v) .* fv(v), v0, vmax, ...
        'RelTol', 1e-8, 'AbsTol', 1e-11);
end
