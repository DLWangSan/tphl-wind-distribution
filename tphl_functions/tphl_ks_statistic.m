function ks = tphl_ks_statistic(data, alpha, mu, lambda, v0)
%TPHL_KS_STATISTIC Kolmogorov-Smirnov test statistic
%   ks = TPHL_KS_STATISTIC(data, alpha, mu, lambda, v0)
%   
%   Returns KS distance between empirical and fitted CDF
%   KS = sup |F_emp(v) - F_fit(v)|
%
%   Reference: Wang et al. (2025), Equation (19)
    
    if nargin < 5, v0 = 0; end
    data_sorted = sort(data);
    n = length(data_sorted);
    p_emp = ((1:n)' - 0.5) / n;
    p_fit = tphl_cdf(data_sorted, alpha, mu, lambda, v0);
    ks = max(abs(p_emp - p_fit));
end

