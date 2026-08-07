%% ENERGY ANALYSIS FIGURE (4-PARAMETER, LITERATURE BOUNDS)
% Uses validated 4-parameter model with k_conc = 4e-6
% Combined RMSE = 0.018986 V, R² = 0.986711

clear; clc; close all;

% ========== MODEL PARAMETERS (4-PARAMETER, LITERATURE BOUNDS) ==========
E_rev = 1.42;
alpha = 0.448938;           % Updated
j0 = 5.428292e-02;          % Updated
R_ohm = 0.169475;           % Updated
k_conc = 4e-6;              % Updated (literature-bounded)
n_e = 2;
R = 8.314;
F = 96485;
T = 348.15;

% ========== MODEL FUNCTION ==========
model = @(j) E_rev + ...
    (R * T / (alpha * n_e * F)) * asinh(j / (2 * j0)) + ...
    R_ohm * j + k_conc * j.^2;

% ========== CURRENT DENSITY RANGE ==========
j_range = linspace(0.01, 2.0, 500);
V_range = model(j_range);

% ========== ENERGY CALCULATIONS ==========
N_cells = 80;
A_cell = 424;  % cm²
I = j_range * A_cell;
P_stack = N_cells * V_range .* I;
mass_H2 = (N_cells * I) / (2 * F) * 0.002016;  % kg/s
specific_energy = (P_stack / 1000) ./ (mass_H2 * 3600);  % kWh/kg
efficiency = 33.3 ./ specific_energy * 100;  % LHV efficiency
energy_per_vol = specific_energy * 0.0899;   % kWh/Nm³

% ========== DISPLAY AT J = 1.0 A/cm² ==========
idx_1A = find(j_range >= 1.0, 1, 'first');
fprintf('\n========================================\n');
fprintf('ENERGY ANALYSIS AT J = 1.0 A/cm²\n');
fprintf('(4-Parameter Model, k_conc = 4e-6)\n');
fprintf('========================================\n');
fprintf('  Specific energy: %.2f kWh/kg\n', specific_energy(idx_1A));
fprintf('  LHV Efficiency: %.2f %%\n', efficiency(idx_1A));
fprintf('  Stack Power: %.2f kW\n', P_stack(idx_1A)/1000);

% ========== CREATE FIGURE ==========
fig = figure('Position', [100, 100, 1000, 700], 'Color', 'white');

% ----- Subplot 1: Power -----
subplot(2, 2, 1);
plot(j_range, P_stack/1000, 'b-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Stack Power [kW]', 'FontSize', 10);
title('(a) Electrical Power Consumption', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);

% ----- Subplot 2: Specific Energy -----
subplot(2, 2, 2);
plot(j_range, specific_energy, 'r-', 'LineWidth', 2.5);
hold on;
yline(33.3, 'k--', 'LHV (33.3 kWh/kg)', 'LineWidth', 1.5);
yline(39.4, 'k:', 'Thermodynamic Min (39.4 kWh/kg)', 'LineWidth', 1.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Specific Energy [kWh/kg H₂]', 'FontSize', 10);
title('(b) Specific Energy Consumption', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 8);
grid on; box on;
xlim([0, 2.0]);

% ----- Subplot 3: Efficiency -----
subplot(2, 2, 3);
plot(j_range, efficiency, 'g-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('LHV Efficiency [%]', 'FontSize', 10);
title('(c) LHV-Based Energy Efficiency', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);

% ----- Subplot 4: Energy per Volume -----
subplot(2, 2, 4);
plot(j_range, energy_per_vol, 'm-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Energy per Volume [kWh/Nm³ H₂]', 'FontSize', 10);
title('(d) Energy Use per Unit Volume', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);

% ----- Main Title -----
sgtitle('Energy Analysis of PEM Electrolyzer (4-Parameter Model)', ...
    'FontSize', 14, 'FontWeight', 'bold');

% ========== SAVE FIGURE ==========
saveas(fig, 'Energy_Analysis_of_PEM_Electrolyzer.jpg');
fprintf('\n✅ Figure saved: Energy_Analysis_of_PEM_Electrolyzer.jpg\n');

% Save as PDF as well
print('Energy_Analysis_of_PEM_Electrolyzer', '-dpdf', '-bestfit');
fprintf('✅ PDF saved: Energy_Analysis_of_PEM_Electrolyzer.pdf\n');
