function [theta, loglik] = tphl_mle(data, seeds, lb, ub, options)
%TPHL_MLE Maximum likelihood estimation for TPHL parameters
%   [theta, loglik] = TPHL_MLE(data, seeds, lb, ub, options)
%   
%   Uses multi-start strategy to find global optimum
%   
%   Inputs:
%     data: wind speed observations (vector)
%     seeds: initial parameter guesses [N x 3] matrix [alpha, mu, lambda]
%            (optional, default seeds will be used if empty)
%     lb: lower bounds [alpha, mu, lambda] (optional)
%     ub: upper bounds [alpha, mu, lambda] (optional)
%     options: optimization options (statset, optional)
%   
%   Outputs:
%     theta: best parameter estimate [alpha, mu, lambda]
%     loglik: log-likelihood value
%
%   Example:
%     v = randn(1000,1)*2 + 3; v = v(v>0);  % Sample wind speeds
%     seeds = [1.0, 2.0, 0.6; 0.7, 2.0, 0.3];
%     lb = [1e-4, -20, 1e-3];
%     ub = [50, 50, 10];
%     opts = statset('MaxIter', 1e4, 'Display', 'off');
%     [theta, LL] = tphl_mle(v, seeds, lb, ub, opts);
%
%   Reference: Wang et al. (2025), Section "Estimation, diagnostics, and 
%   comparative evaluation"
    
    if nargin < 2 || isempty(seeds)
        % Default seeds based on data statistics
        medv = median(data);
        stdv = std(data);
        seeds = [
            1.0, medv, 0.6;
            0.7, medv, 0.3;
            1.5, medv, 1.0;
            0.3, medv, 2.0
        ];
    end
    
    if nargin < 3 || isempty(lb)
        lb = [1e-4, -20, 1e-3];
    end
    
    if nargin < 4 || isempty(ub)
        ub = [50, 50, 10];
    end
    
    if nargin < 5 || isempty(options)
        options = statset('MaxIter', 1e4, 'MaxFunEvals', 2e5, ...
            'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'off');
    end
    
    best_loglik = -inf;
    best_theta = seeds(1, :);
    
    pdf_handle = @(x, a, m, l) tphl_pdf(x, a, m, l, 0);
    
    for i = 1:size(seeds, 1)
        s0 = seeds(i, :);
        try
            theta = mle(data, 'pdf', pdf_handle, ...
                'start', s0, 'LowerBound', lb, 'UpperBound', ub, ...
                'Options', options);
            loglik = sum(log(pdf_handle(data, theta(1), theta(2), theta(3)) + eps));
            if loglik > best_loglik
                best_loglik = loglik;
                best_theta = theta;
            end
        catch
            % Skip failed starting points
        end
    end
    
    theta = best_theta;
    loglik = best_loglik;
end

