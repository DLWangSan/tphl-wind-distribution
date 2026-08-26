function S = phl_survival(v, alpha, mu, lambda)
%PHL_SURVIVAL Survival function of PHL distribution
%   S = PHL_SURVIVAL(v, alpha, mu, lambda)
%   
%   Returns S(v) = 1 - F(v) = (1 + exp(lambda*(v-mu)))^(-alpha)
%
%   Reference: Wang et al. (2025)

    S = (1 + exp(lambda .* (v - mu))).^(-alpha);
end

