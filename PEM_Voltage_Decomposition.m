%% GENERATE VOLTAGE DECOMPOSITION FIGURE (k_conc = 0)
% This script generates a 2-panel voltage decomposition figure with:
%   Subplot 1: Cumulative voltage components (E_rev, η_act, η_ohm, η_conc)
%   Subplot 4: Ohmic and concentration overpotentials
% Subplots 2 and 3 are OMITTED for clarity.

clear; clc; close all;

% ========== MODEL PARAMETERS (k_conc = 0) ==========
E_rev = 1.42;
alpha = 0.449003;
j0 = 5.427630e-02;
R_ohm = 0.169489;
k_conc = 0;
n_e = 2;
R = 8.314;
F = 96485;
T = 348.15;  % 75°C

% ========== MODEL FUNCTION ==========
model = @(j) E_rev + ...
    (R * T / (alpha * n_e * F)) * asinh(j / (2 * j0)) + ...
    R_ohm * j + k_conc * j.^2;

% ========== CURRENT DENSITY RANGE ==========
j_plot = linspace(0.01, 2.0, 200);
V_cell = model(j_plot);

% ========== COMPUTE COMPONENTS ==========
E_rev_T = E_rev * ones(size(j_plot));
eta_act = (R * T / (alpha * n_e * F)) * asinh(j_plot / (2 * j0));
eta_ohm = R_ohm * j_plot;
eta_conc = k_conc * j_plot.^2;  % ZERO

% ========== DISPLAY VERIFICATION ==========
idx_1A = find(j_plot >= 1.0, 1, 'first');
E_rev_1A = E_rev;
eta_act_1A = (R * T / (alpha * n_e * F)) * asinh(j_plot(idx_1A) / (2 * j0));
eta_ohm_1A = R_ohm * j_plot(idx_1A);
eta_conc_1A = 0;
V_cell_1A = E_rev_1A + eta_act_1A + eta_ohm_1A + eta_conc_1A;

fprintf('\n========================================\n');
fprintf('VOLTAGE DECOMPOSITION AT J = 1.0 A/cm²\n');
fprintf('========================================\n');
fprintf('  E_rev     = %.4f V (%.1f%%)\n', E_rev_1A, E_rev_1A/V_cell_1A*100);
fprintf('  η_act     = %.4f V (%.1f%%)\n', eta_act_1A, eta_act_1A/V_cell_1A*100);
fprintf('  η_ohm     = %.4f V (%.1f%%)\n', eta_ohm_1A, eta_ohm_1A/V_cell_1A*100);
fprintf('  η_conc    = %.4f V (%.1f%%)\n', eta_conc_1A, eta_conc_1A/V_cell_1A*100);
fprintf('  V_cell    = %.4f V\n', V_cell_1A);

% ========== CREATE 2-PANEL FIGURE ==========
fig = figure('Position', [100, 100, 1000, 500], 'Color', 'white');

% ----- Subplot 1: Cumulative components (full width, left) -----
subplot(1, 2, 1);
hold on; box on; grid on;

% Stacked area plot
fill([j_plot, fliplr(j_plot)], [E_rev_T, zeros(size(E_rev_T))], ...
    [0.95, 0.95, 0.95], 'EdgeColor', 'none');
fill([j_plot, fliplr(j_plot)], [E_rev_T + eta_act, fliplr(E_rev_T)], ...
    [0.75, 0.85, 0.95], 'EdgeColor', 'none');
fill([j_plot, fliplr(j_plot)], [E_rev_T + eta_act + eta_ohm, fliplr(E_rev_T + eta_act)], ...
    [0.55, 0.70, 0.85], 'EdgeColor', 'none');
fill([j_plot, fliplr(j_plot)], [V_cell, fliplr(E_rev_T + eta_act + eta_ohm)], ...
    [0.35, 0.55, 0.75], 'EdgeColor', 'none');

% Total voltage line
plot(j_plot, V_cell, 'k-', 'LineWidth', 2.5);

xlabel('Current Density [A/cm²]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage [V]', 'FontSize', 12, 'FontWeight', 'bold');
title('(a) Cumulative Voltage Components', 'FontSize', 13, 'FontWeight', 'bold');
legend('E_{rev}', 'η_{act}', 'η_{ohm}', 'η_{conc}', 'V_{cell}', ...
    'Location', 'northwest', 'FontSize', 9);
xlim([0, 2.0]);
ylim([0, 2.5]);

% ----- Subplot 4: Ohmic and Concentration (right) -----
subplot(1, 2, 2);
hold on; box on; grid on;
plot(j_plot, eta_ohm, 'b-', 'LineWidth', 2.5, 'DisplayName', 'η_{ohm} = R_{ohm}·J');
plot(j_plot, eta_conc, 'r-', 'LineWidth', 2.5, 'DisplayName', 'η_{conc} = k_{conc}·J² (0)');
xlabel('Current Density [A/cm²]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Overpotential [V]', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) Ohmic and Concentration Overpotentials', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
xlim([0, 2.0]);
ylim([0, 0.4]);

% ----- Main Title -----
sgtitle('Physics-Based Voltage Decomposition (k_{conc} = 0)', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ----- Save -----
drawnow;
saveas(fig, 'Physics_Based_Voltage_Decomposition_Simplified.jpg');
fprintf('\n✅ Figure saved: Physics_Based_Voltage_Decomposition_Simplified.jpg\n');
