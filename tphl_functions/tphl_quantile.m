function vq = tphl_quantile(q, alpha, mu, lambda, v0)
%TPHL_QUANTILE Quantile function of TPHL distribution
%   vq = TPHL_QUANTILE(q, alpha, mu, lambda, v0)
%   
%   Returns v such that F(v) = q, where q in (0,1)
%
%   Reference: Wang et al. (2025)

    if nargin < 5, v0 = 0; end
    q = min(max(q, 1e-8), 1-1e-8);  % Avoid endpoints
    S0 = phl_survival(v0, alpha, mu, lambda);
    S = max((1 - q) .* S0, eps);
    t = S.^(-1/alpha) - 1;
    vq = mu + (1./lambda) .* log(max(t, eps));
    vq = max(vq, v0);
end

