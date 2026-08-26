function F = tphl_cdf(v, alpha, mu, lambda, v0)
%TPHL_CDF Cumulative distribution function of TPHL distribution
%   F = TPHL_CDF(v, alpha, mu, lambda, v0)
%   
%   Returns F(v) = P(V <= v | V >= v0)
%
%   Reference: Wang et al. (2025), Equation (4)

    if nargin < 5, v0 = 0; end
    S0 = max(phl_survival(v0, alpha, mu, lambda), eps);
    F = zeros(size(v));
    mask = v >= v0;
    S_v = phl_survival(v(mask), alpha, mu, lambda);
    F(mask) = (S0 - S_v) ./ S0;
    F = max(min(F, 1-1e-12), 0);
end

