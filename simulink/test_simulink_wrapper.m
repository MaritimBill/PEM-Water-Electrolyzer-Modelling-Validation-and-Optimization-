%% test_simulink_wrapper.m -- sanity check
clear; clc; close all;

%% ===================== PATH BOOTSTRAP =====================
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
platemoRoot = fileparts(scriptDir);
if exist(fullfile(platemoRoot,'Algorithms'),'dir') ~= 7
    platemoRoot = scriptDir;
end
addpath(genpath(platemoRoot));
fprintf('Path root: %s\n\n', platemoRoot);

%% ===================== MODEL =====================
if exist('PEM_Plant.slx', 'file') ~= 4
    fprintf('PEM_Plant.slx missing -> rebuilding.\n');
    build_PEM_Simulink();
end

w = PEM_Simulink_Wrapper();
P = PEM_Problem_Scaled();
P.Setting();

x = [100, 20, 10000, 70, 10e5];

metrics = P.CalculateMetrics(x);
[V_sl, P_sl, m_sl] = w.simulate(x(1), x(2), x(3), x(4), x(5));

%% ===================== REPORT =====================
fprintf('\n=== Analytical vs Simulink ===\n');
fprintf('V_cell  :  ana=%.4f V   sim=%.4f V   rel=%.4f%%\n', ...
    metrics.V_cell, V_sl, 100*abs(metrics.V_cell - V_sl)/max(metrics.V_cell,1e-12));
fprintf('P_stack :  ana=%.3f W   sim=%.3f W   rel=%.4f%%\n', ...
    metrics.P_stack, P_sl, 100*abs(metrics.P_stack - P_sl)/max(metrics.P_stack,1));
fprintf('m_H2    :  ana=%.6e kg/s  sim=%.6e kg/s  rel=%.4f%%\n', ...
    metrics.mass_H2, m_sl, 100*abs(metrics.mass_H2 - m_sl)/max(metrics.mass_H2,1e-12));

fprintf('\nTiming (fast-restart warm):\n');
t = tic;
for k = 1:20
    w.simulate(x(1), x(2), x(3)+k, x(4), x(5));
end
dt = toc(t);
fprintf('  20 sims in %.3f s -> %.2f ms/sim\n', dt, dt/20*1000);

delete(w);
fprintf('\nOK.\n');
