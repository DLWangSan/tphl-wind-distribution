function f = phl_pdf(v, alpha, mu, lambda)
%PHL_PDF Probability density function of PHL distribution
%   f = PHL_PDF(v, alpha, mu, lambda)
%   
%   Inputs:
%     v: wind speed (m/s)
%     alpha: shape parameter (>0)
%     mu: location parameter
%     lambda: scale parameter (>0)
%   
%   Output:
%     f: PDF value
%
%   Reference: Wang et al. (2026) "A Physically Consistent TPHL Distribution 
%   for Wind-Speed Modeling in Coastal-Inland Regimes" Ocean-Land-Atmosphere Research, doi:10.34133/olar.0163

    z = lambda .* (v - mu);
    f = alpha .* lambda .* exp(z) ./ (1 + exp(z)).^(alpha + 1);
end

