%% test_simulink_wrapper.m -- sanity check, run me first
clear; clc;

if exist('PEM_Plant.slx', 'file') ~= 4
    build_PEM_Simulink();
end

w = PEM_Simulink_Wrapper();

% Reference analytical evaluation via your own PEM_Problem_Scaled
P = PEM_Problem_Scaled();
P.Setting();

% Test vector well inside feasible region
x = [100, 20, 10000, 70, 10e5];   % cm2, -, A/m2, C, Pa

metrics = P.CalculateMetrics(x);
[V_sl, P_sl, m_sl, T_sl] = w.simulate(x(1), x(2), x(3), x(4), x(5));

fprintf('\n=== Analytical vs Simulink ===\n');
fprintf('V_cell  :  ana=%.4f V   sim=%.4f V   rel=%.4f%%\n', ...
    metrics.V_cell, V_sl, 100*abs(metrics.V_cell - V_sl)/metrics.V_cell);
fprintf('P_stack :  ana=%.3f W   sim=%.3f W   rel=%.4f%%\n', ...
    metrics.P_stack, P_sl, 100*abs(metrics.P_stack - P_sl)/max(metrics.P_stack,1));
fprintf('m_H2    :  ana=%.6e kg/s  sim=%.6e kg/s  rel=%.4f%%\n', ...
    metrics.mass_H2, m_sl, 100*abs(metrics.mass_H2 - m_sl)/max(metrics.mass_H2,1e-12));
fprintf('T_cell  :  sim=%.2f K  (analytical model has no thermal state)\n', T_sl);

fprintf('\nTiming (fast-restart warm):\n');
t = tic;
for k = 1:20
    w.simulate(x(1), x(2), x(3)+k, x(4), x(5));
end
fprintf('  20 sims in %.3f s -> %.2f ms/sim\n', toc(t), toc(t)/20*1000);

delete(w);
fprintf('\nOK.\n');