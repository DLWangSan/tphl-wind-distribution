function aic = tphl_aic(data, alpha, mu, lambda, v0)
%TPHL_AIC Akaike Information Criterion
%   aic = TPHL_AIC(data, alpha, mu, lambda, v0)
%   
%   AIC = 2k - 2*log(L), where k=3 for TPHL
%
%   Reference: Wang et al. (2025), Equation (21)
    
    if nargin < 5, v0 = 0; end
    loglik = sum(log(tphl_pdf(data, alpha, mu, lambda, v0) + eps));
    aic = 2*3 - 2*loglik;
end

