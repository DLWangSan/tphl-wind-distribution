function [meanP, varP, skewP, kurtP] = tphl_power_stats(alpha, mu, lambda, c, v0)
%TPHL_POWER_STATS Compute wind power statistics from TPHL parameters
%   [meanP, varP, skewP, kurtP] = TPHL_POWER_STATS(alpha, mu, lambda, c, v0)
%   
%   Computes statistics of P = c * V^3 from TPHL parameters
%   
%   Inputs:
%     alpha, mu, lambda: TPHL parameters
%     c: power constant (default: 1.0, or use 0.5*rho for physical power)
%     v0: truncation point (default: 0)
%   
%   Outputs:
%     meanP: mean power
%     varP: variance of power
%     skewP: skewness of power
%     kurtP: excess kurtosis of power
%
%   Reference: Wang et al. (2025), Equations (13)-(16)
    
    if nargin < 4, c = 1.0; end
    if nargin < 5, v0 = 0; end
    
    % Determine integration upper bound
    v99 = tphl_quantile(0.999, alpha, mu, lambda, v0);
    vmax = max(v99, v0) + 10 / max(lambda, 1e-3);
    
    fv = @(v) tphl_pdf(v, alpha, mu, lambda, v0);
    
    % Compute moments E[V^m] for m = 3, 6, 9, 12
    Ev3  = integral(@(v) v.^3  .* fv(v), v0, vmax, 'RelTol', 1e-7, 'AbsTol', 1e-10);
    Ev6  = integral(@(v) v.^6  .* fv(v), v0, vmax, 'RelTol', 1e-7, 'AbsTol', 1e-10);
    Ev9  = integral(@(v) v.^9  .* fv(v), v0, vmax, 'RelTol', 1e-6, 'AbsTol', 1e-9);
    Ev12 = integral(@(v) v.^12 .* fv(v), v0, vmax, 'RelTol', 1e-6, 'AbsTol', 1e-9);
    
    % Power moments: P = c * V^3
    m1 = c * Ev3;
    m2 = c^2 * Ev6;
    m3 = c^3 * Ev9;
    m4 = c^4 * Ev12;
    
    % Compute statistics
    meanP = m1;
    varP = max(m2 - m1^2, 0);
    sdP = sqrt(max(varP, 0));
    
    if sdP > 0
        skewP = (m3 - 3*m1*m2 + 2*m1^3) / sdP^3;
        kurtP = (m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4) / sdP^4 - 3;
    else
        skewP = NaN;
        kurtP = NaN;
    end
end

