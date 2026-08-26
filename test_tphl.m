%% Quick Test Script for TPHL Functions
% This script performs a quick test to ensure all functions work correctly

clear; clc;

fprintf('=== Testing TPHL Functions ===\n\n');

% Add function path
addpath('tphl_functions');

% Test parameters
alpha = 0.3;
mu = 2.0;
lambda = 1.5;
v_test = [0.5, 1.0, 2.0, 3.0, 5.0];

%% Test 1: PDF
fprintf('Test 1: PDF function... ');
try
    pdf_vals = tphl_pdf(v_test, alpha, mu, lambda);
    if all(isfinite(pdf_vals)) && all(pdf_vals >= 0)
        fprintf('PASS\n');
    else
        fprintf('FAIL (invalid values)\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
end

%% Test 2: CDF
fprintf('Test 2: CDF function... ');
try
    cdf_vals = tphl_cdf(v_test, alpha, mu, lambda);
    if all(isfinite(cdf_vals)) && all(cdf_vals >= 0) && all(cdf_vals <= 1)
        fprintf('PASS\n');
    else
        fprintf('FAIL (invalid values)\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
end

%% Test 3: Quantile
fprintf('Test 3: Quantile function... ');
try
    q_test = [0.1, 0.5, 0.9];
    vq = tphl_quantile(q_test, alpha, mu, lambda);
    if all(isfinite(vq)) && all(vq >= 0)
        fprintf('PASS\n');
    else
        fprintf('FAIL (invalid values)\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
end

%% Test 4: MLE (with synthetic data)
fprintf('Test 4: MLE function... ');
try
    % Generate synthetic data
    n = 500;
    u = rand(n, 1);
    S0 = (1 + exp(-lambda * mu))^(-alpha);
    S = (1 - u) .* S0;
    t = S.^(-1/alpha) - 1;
    v_synth = mu + (1/lambda) * log(max(t, eps));
    v_synth = max(v_synth, 0);
    
    seeds = [1.0, median(v_synth), 0.6; 0.7, median(v_synth), 0.3];
    lb = [1e-4, -20, 1e-3];
    ub = [50, 50, 10];
    opts = statset('MaxIter', 1000, 'Display', 'off');
    
    [theta, loglik] = tphl_mle(v_synth, seeds, lb, ub, opts);
    if all(isfinite(theta)) && isfinite(loglik)
        fprintf('PASS (alpha=%.3f, mu=%.3f, lambda=%.3f)\n', theta(1), theta(2), theta(3));
    else
        fprintf('FAIL (invalid values)\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
end

%% Test 5: Power Statistics
fprintf('Test 5: Power statistics... ');
try
    [meanP, varP, skewP, kurtP] = tphl_power_stats(alpha, mu, lambda, 1.0);
    if all(isfinite([meanP, varP, skewP, kurtP]))
        fprintf('PASS (mean=%.3f, var=%.3f, skew=%.3f, kurt=%.3f)\n', ...
            meanP, varP, skewP, kurtP);
    else
        fprintf('FAIL (invalid values)\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
end

%% Test 6: Goodness-of-fit metrics
fprintf('Test 6: Goodness-of-fit metrics... ');
try
    % Use synthetic data from Test 4
    if exist('v_synth', 'var')
        ks = tphl_ks_statistic(v_synth, alpha, mu, lambda);
        aic = tphl_aic(v_synth, alpha, mu, lambda);
        bic = tphl_bic(v_synth, alpha, mu, lambda);
        if all(isfinite([ks, aic, bic]))
            fprintf('PASS (KS=%.4f, AIC=%.2f, BIC=%.2f)\n', ks, aic, bic);
        else
            fprintf('FAIL (invalid values)\n');
        end
    else
        fprintf('SKIP (no test data)\n');
    end
catch ME
    fprintf('FAIL: %s\n', ME.message);
end

fprintf('\n=== All tests completed ===\n');

