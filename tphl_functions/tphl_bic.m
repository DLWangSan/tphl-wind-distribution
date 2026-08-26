function bic = tphl_bic(data, alpha, mu, lambda, v0)
%TPHL_BIC Bayesian Information Criterion
%   bic = TPHL_BIC(data, alpha, mu, lambda, v0)
%   
%   BIC = k*log(n) - 2*log(L), where k=3 for TPHL
%
%   Reference: Wang et al. (2025), Equation (21)
    
    if nargin < 5, v0 = 0; end
    n = length(data);
    loglik = sum(log(tphl_pdf(data, alpha, mu, lambda, v0) + eps));
    bic = 3*log(n) - 2*loglik;
end

