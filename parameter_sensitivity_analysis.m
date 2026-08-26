%% Parameter Sensitivity Analysis for TPHL Distribution
% This script analyzes how TPHL parameters affect the distribution shape
% and wind power statistics.
%
% Reference: Wang et al. (2026) "A Physically Consistent TPHL Distribution 
% for Wind-Speed Modeling in Coastal-Inland Regimes" Ocean-Land-Atmosphere Research, doi:10.34133/olar.0163

clear; clc; close all;

% Add function path
addpath('tphl_functions');

% Create figures directory if it doesn't exist
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% 1. Base parameters (typical inland station)
alpha_base = 0.3;
mu_base = 1.8;
lambda_base = 4.0;

fprintf('=== Parameter Sensitivity Analysis ===\n');
fprintf('Base parameters:\n');
fprintf('  alpha  = %.2f\n', alpha_base);
fprintf('  mu     = %.2f m/s\n', mu_base);
fprintf('  lambda = %.2f\n', lambda_base);

%% 2. Sensitivity to alpha (tail heaviness)
fprintf('\n=== Sensitivity to alpha (tail heaviness) ===\n');
alpha_range = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0];
v_plot = linspace(0, 10, 500);

figure('Position', [100, 100, 1400, 500]);

subplot(1, 3, 1);
hold on;
colors = lines(length(alpha_range));
for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    pdf_vals = tphl_pdf(v_plot, alpha, mu_base, lambda_base);
    plot(v_plot, pdf_vals, '-', 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('α=%.2f', alpha));
end
xlabel('Wind speed (m/s)', 'FontSize', 12);
ylabel('PDF', 'FontSize', 12);
title('Effect of α on PDF', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

% Power statistics vs alpha
subplot(1, 3, 2);
meanP_alpha = zeros(size(alpha_range));
skewP_alpha = zeros(size(alpha_range));
for i = 1:length(alpha_range)
    [meanP, ~, skewP, ~] = tphl_power_stats(alpha_range(i), mu_base, lambda_base, 1.0);
    meanP_alpha(i) = meanP;
    skewP_alpha(i) = skewP;
end
yyaxis left;
plot(alpha_range, meanP_alpha, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Mean power', 'FontSize', 12);
yyaxis right;
plot(alpha_range, skewP_alpha, 's-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Power skewness', 'FontSize', 12);
xlabel('α (shape parameter)', 'FontSize', 12);
title('Power Statistics vs α', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('Mean power', 'Power skewness', 'Location', 'best');

% Tail behavior
subplot(1, 3, 3);
p_tail = [0.95, 0.99, 0.999];
v_tail = zeros(length(alpha_range), length(p_tail));
for i = 1:length(alpha_range)
    for j = 1:length(p_tail)
        v_tail(i, j) = tphl_quantile(p_tail(j), alpha_range(i), mu_base, lambda_base);
    end
end
plot(alpha_range, v_tail(:,1), 'o-', 'LineWidth', 2, 'DisplayName', 'P95');
hold on;
plot(alpha_range, v_tail(:,2), 's-', 'LineWidth', 2, 'DisplayName', 'P99');
plot(alpha_range, v_tail(:,3), '^-', 'LineWidth', 2, 'DisplayName', 'P99.9');
xlabel('α (shape parameter)', 'FontSize', 12);
ylabel('Quantile (m/s)', 'FontSize', 12);
title('Tail Quantiles vs α', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

sgtitle('Sensitivity to α (Tail Heaviness)', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'figures/sensitivity_alpha.png', 'png');
saveas(gcf, 'figures/sensitivity_alpha.fig', 'fig');
fprintf('Figure saved: figures/sensitivity_alpha.png\n');

%% 3. Sensitivity to mu (location parameter)
fprintf('\n=== Sensitivity to mu (location parameter) ===\n');
mu_range = [0.5, 1.0, 1.8, 2.5, 3.5, 5.0];

figure('Position', [100, 100, 1400, 500]);

subplot(1, 3, 1);
hold on;
colors = lines(length(mu_range));
for i = 1:length(mu_range)
    mu = mu_range(i);
    pdf_vals = tphl_pdf(v_plot, alpha_base, mu, lambda_base);
    plot(v_plot, pdf_vals, '-', 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('μ=%.2f', mu));
end
xlabel('Wind speed (m/s)', 'FontSize', 12);
ylabel('PDF', 'FontSize', 12);
title('Effect of μ on PDF', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

% Mean wind speed vs mu
subplot(1, 3, 2);
meanV_mu = zeros(size(mu_range));
for i = 1:length(mu_range)
    % Compute mean wind speed numerically
    v99 = tphl_quantile(0.999, alpha_base, mu_range(i), lambda_base);
    vmax = max(v99, 0) + 10 / max(lambda_base, 1e-3);
    fv = @(v) tphl_pdf(v, alpha_base, mu_range(i), lambda_base, 0);
    meanV_mu(i) = integral(@(v) v .* fv(v), 0, vmax, 'RelTol', 1e-7);
end
plot(mu_range, meanV_mu, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('μ (location parameter)', 'FontSize', 12);
ylabel('Mean wind speed (m/s)', 'FontSize', 12);
title('Mean Wind Speed vs μ', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Power statistics vs mu
subplot(1, 3, 3);
meanP_mu = zeros(size(mu_range));
for i = 1:length(mu_range)
    [meanP, ~, ~, ~] = tphl_power_stats(alpha_base, mu_range(i), lambda_base, 1.0);
    meanP_mu(i) = meanP;
end
plot(mu_range, meanP_mu, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('μ (location parameter)', 'FontSize', 12);
ylabel('Mean power', 'FontSize', 12);
title('Mean Power vs μ', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

sgtitle('Sensitivity to μ (Location Parameter)', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'figures/sensitivity_mu.png', 'png');
saveas(gcf, 'figures/sensitivity_mu.fig', 'fig');
fprintf('Figure saved: figures/sensitivity_mu.png\n');

%% 4. Sensitivity to lambda (scale parameter)
fprintf('\n=== Sensitivity to lambda (scale parameter) ===\n');
lambda_range = [1.0, 2.0, 4.0, 6.0, 8.0];

figure('Position', [100, 100, 1400, 500]);

subplot(1, 3, 1);
hold on;
colors = lines(length(lambda_range));
for i = 1:length(lambda_range)
    lambda = lambda_range(i);
    pdf_vals = tphl_pdf(v_plot, alpha_base, mu_base, lambda);
    plot(v_plot, pdf_vals, '-', 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('λ=%.2f', lambda));
end
xlabel('Wind speed (m/s)', 'FontSize', 12);
ylabel('PDF', 'FontSize', 12);
title('Effect of λ on PDF', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

% Distribution width vs lambda
subplot(1, 3, 2);
stdV_lambda = zeros(size(lambda_range));
for i = 1:length(lambda_range)
    % Compute std wind speed numerically
    v99 = tphl_quantile(0.999, alpha_base, mu_base, lambda_range(i));
    vmax = max(v99, 0) + 10 / max(lambda_range(i), 1e-3);
    fv = @(v) tphl_pdf(v, alpha_base, mu_base, lambda_range(i), 0);
    meanV = integral(@(v) v .* fv(v), 0, vmax, 'RelTol', 1e-7);
    meanV2 = integral(@(v) v.^2 .* fv(v), 0, vmax, 'RelTol', 1e-7);
    stdV_lambda(i) = sqrt(max(meanV2 - meanV^2, 0));
end
plot(lambda_range, stdV_lambda, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('λ (scale parameter)', 'FontSize', 12);
ylabel('Std wind speed (m/s)', 'FontSize', 12);
title('Std Wind Speed vs λ', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Power variance vs lambda
subplot(1, 3, 3);
varP_lambda = zeros(size(lambda_range));
for i = 1:length(lambda_range)
    [~, varP, ~, ~] = tphl_power_stats(alpha_base, mu_base, lambda_range(i), 1.0);
    varP_lambda(i) = varP;
end
plot(lambda_range, varP_lambda, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('λ (scale parameter)', 'FontSize', 12);
ylabel('Power variance', 'FontSize', 12);
title('Power Variance vs λ', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

sgtitle('Sensitivity to λ (Scale Parameter)', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'figures/sensitivity_lambda.png', 'png');
saveas(gcf, 'figures/sensitivity_lambda.fig', 'fig');
fprintf('Figure saved: figures/sensitivity_lambda.png\n');

%% 5. Summary table
fprintf('\n=== Summary: Parameter Effects ===\n');
fprintf('Parameter | Effect on Distribution | Effect on Power\n');
fprintf('----------|----------------------|----------------\n');
fprintf('α (alpha) | Tail heaviness       | Skewness, extremes\n');
fprintf('μ (mu)    | Central location     | Mean power\n');
fprintf('λ (lambda)| Distribution width   | Variance\n');

fprintf('\nAnalysis complete!\n');

