%% Example: TPHL Distribution Fitting for Wind Speed Data
% This script demonstrates how to fit the TPHL distribution to wind speed data
% and compare it with the Weibull distribution.
%
% Reference: Wang et al. (2026) "A Physically Consistent TPHL Distribution 
% for Wind-Speed Modeling in Coastal-Inland Regimes" Ocean-Land-Atmosphere Research, doi:10.34133/olar.0163

clear; clc; close all;

% Add function path
addpath('tphl_functions');

%% 1. Generate or load sample wind speed data
% Option A: Generate synthetic data from TPHL
alpha_true = 0.3;
mu_true = 2.5;
lambda_true = 1.5;
n_samples = 1000;

% Generate random samples (using inverse CDF method)
u = rand(n_samples, 1);
v0 = 0;
S0 = (1 + exp(-lambda_true * mu_true))^(-alpha_true);
S = (1 - u) .* S0;
t = S.^(-1/alpha_true) - 1;
v_synthetic = mu_true + (1/lambda_true) * log(max(t, eps));
v_synthetic = max(v_synthetic, v0);

% Option B: Load actual wind speed data
% load('wind_data.mat', 'v');  % Load from file
% v = v(v > 0.1);  % Remove very small values

% Use synthetic data for this example
v = v_synthetic;
v = v(v > 0.1);  % Remove very small values

fprintf('Sample size: %d\n', length(v));
fprintf('Mean: %.3f m/s, Std: %.3f m/s\n', mean(v), std(v));

%% 2. Fit TPHL distribution
fprintf('\n=== Fitting TPHL distribution ===\n');

% Set up MLE options
seeds = [
    1.0, median(v), 0.6;
    0.7, median(v), 0.3;
    1.5, median(v), 1.0;
    0.3, median(v), 2.0
];
lb = [1e-4, -20, 1e-3];
ub = [50, 50, 10];
opts = statset('MaxIter', 1e4, 'MaxFunEvals', 2e5, ...
    'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'off');

% Fit TPHL
[theta_tphl, loglik_tphl] = tphl_mle(v, seeds, lb, ub, opts);
alpha_tphl = theta_tphl(1);
mu_tphl = theta_tphl(2);
lambda_tphl = theta_tphl(3);

fprintf('TPHL parameters:\n');
fprintf('  alpha = %.4f\n', alpha_tphl);
fprintf('  mu    = %.4f m/s\n', mu_tphl);
fprintf('  lambda = %.4f\n', lambda_tphl);
fprintf('  Log-likelihood = %.2f\n', loglik_tphl);

%% 3. Fit Weibull distribution (for comparison)
fprintf('\n=== Fitting Weibull distribution ===\n');
try
    pd_wb = fitdist(v, 'Weibull');
    scale_wb = pd_wb.A;
    shape_wb = pd_wb.B;
    loglik_wb = sum(log(pdf(pd_wb, v) + eps));
    fprintf('Weibull parameters:\n');
    fprintf('  scale = %.4f m/s\n', scale_wb);
    fprintf('  shape = %.4f\n', shape_wb);
    fprintf('  Log-likelihood = %.2f\n', loglik_wb);
catch
    warning('Weibull fitting failed');
    scale_wb = NaN;
    shape_wb = NaN;
    loglik_wb = NaN;
end

%% 4. Compute goodness-of-fit metrics
fprintf('\n=== Goodness-of-fit metrics ===\n');

% KS statistic
ks_tphl = tphl_ks_statistic(v, alpha_tphl, mu_tphl, lambda_tphl);
if ~isnan(scale_wb)
    v_sort = sort(v);
    n = length(v_sort);
    p_emp = ((1:n)' - 0.5) / n;
    p_wb = wblcdf(v_sort, scale_wb, shape_wb);
    ks_wb = max(abs(p_emp - p_wb));
else
    ks_wb = NaN;
end

fprintf('KS distance:\n');
fprintf('  TPHL:    %.4f\n', ks_tphl);
fprintf('  Weibull: %.4f\n', ks_wb);

% AIC
aic_tphl = tphl_aic(v, alpha_tphl, mu_tphl, lambda_tphl);
if ~isnan(loglik_wb)
    aic_wb = 2*2 - 2*loglik_wb;  % Weibull has 2 parameters
else
    aic_wb = NaN;
end

fprintf('AIC:\n');
fprintf('  TPHL:    %.2f\n', aic_tphl);
fprintf('  Weibull: %.2f\n', aic_wb);

% BIC
bic_tphl = tphl_bic(v, alpha_tphl, mu_tphl, lambda_tphl);
if ~isnan(loglik_wb)
    bic_wb = 2*log(length(v)) - 2*loglik_wb;
else
    bic_wb = NaN;
end

fprintf('BIC:\n');
fprintf('  TPHL:    %.2f\n', bic_tphl);
fprintf('  Weibull: %.2f\n', bic_wb);

%% 5. Compute wind power statistics
fprintf('\n=== Wind power statistics ===\n');
power_const = 1.0;  % Use 1.0 for normalized power, or 0.5*rho for physical
[meanP, varP, skewP, kurtP] = tphl_power_stats(alpha_tphl, mu_tphl, lambda_tphl, power_const);

fprintf('From TPHL parameters:\n');
fprintf('  Mean power:  %.4f\n', meanP);
fprintf('  Variance:    %.4f\n', varP);
fprintf('  Skewness:    %.4f\n', skewP);
fprintf('  Excess Kurtosis: %.4f\n', kurtP);

% Compare with empirical statistics
P_emp = v.^3;
fprintf('\nEmpirical (from data):\n');
fprintf('  Mean power:  %.4f\n', mean(P_emp));
fprintf('  Variance:    %.4f\n', var(P_emp));
fprintf('  Skewness:    %.4f\n', skewness(P_emp));
fprintf('  Excess Kurtosis: %.4f\n', kurtosis(P_emp) - 3);

%% 6. Plot PDF comparison
fprintf('\n=== Generating plots ===\n');

figure('Position', [100, 100, 1200, 400]);

% PDF plot
subplot(1, 3, 1);
histogram(v, 30, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7], ...
    'EdgeColor', 'none', 'DisplayName', 'Observed');
hold on;
x_plot = linspace(0, max(v)*1.2, 500);
y_tphl = tphl_pdf(x_plot, alpha_tphl, mu_tphl, lambda_tphl);
plot(x_plot, y_tphl, 'r-', 'LineWidth', 2, 'DisplayName', 'TPHL');
if ~isnan(scale_wb)
    y_wb = wblpdf(x_plot, scale_wb, shape_wb);
    plot(x_plot, y_wb, 'b--', 'LineWidth', 2, 'DisplayName', 'Weibull');
end
xlabel('Wind speed (m/s)', 'FontSize', 12);
ylabel('PDF', 'FontSize', 12);
title('PDF Comparison', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

% P-P plot
subplot(1, 3, 2);
v_sort = sort(v);
n = length(v_sort);
p_emp = ((1:n)' - 0.5) / n;
p_tphl = tphl_cdf(v_sort, alpha_tphl, mu_tphl, lambda_tphl);
plot([0 1], [0 1], 'k:', 'LineWidth', 1);
hold on;
plot(p_emp, p_tphl, 'r-', 'LineWidth', 1.5, 'DisplayName', 'TPHL');
if ~isnan(scale_wb)
    p_wb = wblcdf(v_sort, scale_wb, shape_wb);
    plot(p_emp, p_wb, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Weibull');
end
xlabel('Empirical probability', 'FontSize', 12);
ylabel('Model probability', 'FontSize', 12);
title('P-P Plot', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
axis equal;

% Q-Q plot
subplot(1, 3, 3);
p = linspace(0.01, 0.99, 200)';
v_emp = quantile(v, p);
v_tphl_q = tphl_quantile(p, alpha_tphl, mu_tphl, lambda_tphl);
plot([0 max(v_emp)], [0 max(v_emp)], 'k:', 'LineWidth', 1);
hold on;
plot(v_emp, v_tphl_q, 'r-', 'LineWidth', 1.5, 'DisplayName', 'TPHL');
if ~isnan(scale_wb)
    v_wb_q = wblinv(p, scale_wb, shape_wb);
    plot(v_emp, v_wb_q, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Weibull');
end
xlabel('Empirical quantile (m/s)', 'FontSize', 12);
ylabel('Model quantile (m/s)', 'FontSize', 12);
title('Q-Q Plot', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
axis equal;

sgtitle('TPHL vs Weibull Distribution Comparison', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('Plots generated successfully!\n');

