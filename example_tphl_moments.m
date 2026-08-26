%% TPHL 各阶矩与矩母函数（MGF）示例
%
% 说明（与老师强调的“亮点”对应）：
%   - 补充材料 S1 式 (S8): E[V^r] = 截断 PHL 上的原矩
%   - 补充材料 S1 式 (S9): M_V(t) = E[e^{tV}]  —— 矩母函数
%   - 在 t=0 处：M_V^{(k)}(0) = E[V^k]（收敛时），故高阶矩不必每次写复杂积分，
%     也可从 MGF 求导得到；论文同时给出 logistic→Beta 的闭式（S10）。
%
% 本脚本用数值积分实现 S8/S9，并对比“直接算原矩”与“MGF 数值求导”。
%
% Reference: Wang et al. (2026) OLAR, doi:10.34133/olar.0163

clear; clc;
addpath('tphl_functions');

%% 示例参数（可改成你拟合得到的 alpha, mu, lambda）
alpha = 0.35;
mu    = 2.0;
lambda = 1.2;
c_wpd = 0.5 * 1.225;   % WPD 用 c = rho/2

fprintf('=== TPHL moments & MGF demo ===\n');
fprintf('alpha=%.4f  mu=%.4f  lambda=%.4f\n\n', alpha, mu, lambda);

R = tphl_moments_report(alpha, mu, lambda, ...
    'max_r', 6, ...
    'c_power', c_wpd, ...
    'mgf_check', true, ...
    't_grid', linspace(-0.06, 0.06, 25));

fprintf('--- Raw moments E[V^r] (Supp. S8) ---\n');
for r = 1:numel(R.raw.M)
    fprintf('  M%d = E[V^%d] = %.6f\n', r, r, R.raw.M(r));
end

fprintf('\n--- Wind speed shape (Main text Eq. 9-12) ---\n');
fprintf('  mean(V)  = %.6f m/s\n', R.speed.mean);
fprintf('  var(V)   = %.6f (m/s)^2\n', R.speed.var);
fprintf('  skew(V)  = %.6f\n', R.speed.skew);
fprintf('  kurt(V)  = %.6f  (Pearson; normal baseline ~ 3)\n', R.speed.kurt);

fprintf('\n--- Wind power P = c*V^3 (Main text Eq. 13-16, c=rho/2) ---\n');
fprintf('  E[P]     = %.4f W/m^2\n', R.power.meanP);
fprintf('  Var(P)   = %.4f\n', R.power.varP);
fprintf('  skew(P)  = %.4f\n', R.power.skewP);
fprintf('  kurt(P)  = %.4f  (Pearson; paper Eq. 16)\n', R.power.kurtP);

if ~isempty(fieldnames(R.mgf)) && isfield(R.mgf, 'order')
    fprintf('\n--- MGF check: E[V^k] vs d^k M/dt^k|_{t=0} (Supp. S9) ---\n');
    T = table(R.mgf.order, R.mgf.from_raw_moment, R.mgf.from_mgf_derivative, R.mgf.abs_diff, ...
        'VariableNames', {'k', 'E_Vk_direct', 'E_Vk_from_MGF', 'abs_diff'});
    disp(T);
    fprintf('  (k=1,2: direct raw moments vs numerical MGF derivative at t=0; higher k use tphl_raw_moment.)\n');
end

%% 可选：从 HadISD 拟合后再算矩（需本地 nc 路径）
% Uncomment and set path if available:
% v = ... load windspeeds ...
% [theta, ~] = tphl_mle(v, seeds, lb, ub, opts);
% R_site = tphl_moments_report(theta(1), theta(2), theta(3), 'max_r', 12, 'mgf_check', true);

fprintf('\nDone.\n');
