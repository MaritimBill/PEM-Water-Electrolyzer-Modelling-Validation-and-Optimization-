%% test_simulink_wrapper.m
%  Build (if needed), wrap, simulate, compare against analytical, time it.

clc;
fprintf('\n=== test_simulink_wrapper ===\n');

% ---- Ensure the model exists ----
if exist('PEM_Plant.slx','file') ~= 2
    fprintf('PEM_Plant.slx not found, building...\n');
    build_PEM_Simulink();
end

% ---- Test point ----
A = 100;  N = 20;  J = 1.0;  T_C = 70;  P_bar = 10;

% ---- Analytical reference ----
F = 96485;  R = 8.314;  M_H2 = 2.016e-3;
E_rev  = 1.42;   alpha = 0.448938;   j0 = 5.428292e-2;
R_ohm  = 0.169475;  k_conc = 4e-6;   n_e = 2;

T_K      = T_C + 273.15;
E_rev_T  = E_rev - 0.0009*(T_C - 75);
f        = R*T_K/(alpha*n_e*F);
eta_act  = f * asinh(J/(2*j0));
eta_ohm  = R_ohm * J;
eta_conc = k_conc * J^2;
V_nernst = (R*T_K/(2*F)) * log(P_bar);
V_ana    = E_rev_T + eta_act + eta_ohm + eta_conc + V_nernst;

I_stack  = J * A;
P_ana    = N * V_ana * I_stack;
m_ana    = N * I_stack / (2*F) * M_H2;

% ---- Simulate ----
w = PEM_Simulink_Wrapper('PEM_Plant');
[V, P, m] = w.simulate(A, N, J, T_C, P_bar);

% ---- Report ----
fprintf('\n=== Analytical vs Simulink ===\n');
fprintf('%-8s  ana=%-14.6g  sim=%-14.6g  rel=%.4f%%\n', ...
        'V_cell',  V_ana, V, 100*abs(V - V_ana)/abs(V_ana));
fprintf('%-8s  ana=%-14.6g  sim=%-14.6g  rel=%.4f%%\n', ...
        'P_stack', P_ana, P, 100*abs(P - P_ana)/abs(P_ana));
fprintf('%-8s  ana=%-14.6g  sim=%-14.6g  rel=%.4f%%\n', ...
        'm_H2',    m_ana, m, 100*abs(m - m_ana)/abs(m_ana));

% ---- Timing (warm) ----
w.simulate(A, N, J, T_C, P_bar);          % warm
t0 = tic;
n  = 20;
for k = 1:n
    w.simulate(A, N, J, T_C, P_bar);
end
elapsed = toc(t0);
fprintf('\nTiming (accelerator, warm):\n');
fprintf('  %d sims in %.3f s -> %.2f ms/sim\n', ...
        n, elapsed, 1000*elapsed/n);

fprintf('\nOK.\n\n');
