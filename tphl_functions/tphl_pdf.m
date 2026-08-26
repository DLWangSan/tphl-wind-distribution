function f = tphl_pdf(v, alpha, mu, lambda, v0)
%TPHL_PDF Probability density function of TPHL distribution
%   f = TPHL_PDF(v, alpha, mu, lambda, v0)
%   
%   TPHL is PHL truncated at v >= v0 (typically v0=0 for wind speed)
%   
%   Inputs:
%     v: wind speed (m/s), must be >= v0
%     alpha: shape parameter (>0), controls tail heaviness
%     mu: location parameter, characteristic wind speed
%     lambda: scale parameter (>0), controls distribution width
%     v0: truncation point (default: 0)
%   
%   Output:
%     f: PDF value
%
%   Reference: Wang et al. (2025), Equation (6)

    if nargin < 5, v0 = 0; end
    S0 = max(phl_survival(v0, alpha, mu, lambda), eps);
    f = zeros(size(v));
    mask = v >= v0;
    f(mask) = phl_pdf(v(mask), alpha, mu, lambda) ./ S0;
end

