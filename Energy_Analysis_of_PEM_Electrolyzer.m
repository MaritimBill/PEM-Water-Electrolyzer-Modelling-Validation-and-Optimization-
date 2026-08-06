%% ENERGY ANALYSIS FIGURE (k_conc = 0)
clear; clc; close all;

% Parameters
E_rev = 1.42;
alpha = 0.449003;
j0 = 5.427630e-02;
R_ohm = 0.169489;
k_conc = 0;
n_e = 2;
R = 8.314;
F = 96485;
T = 348.15;

% Model
model = @(j) E_rev + ...
    (R * T / (alpha * n_e * F)) * asinh(j / (2 * j0)) + ...
    R_ohm * j + k_conc * j.^2;

% Current density range
j_range = linspace(0.01, 2.0, 500);
V_range = model(j_range);

% Energy calculations
N_cells = 80;
A_cell = 424;  % cm²
I = j_range * A_cell;
P_stack = N_cells * V_range .* I;
mass_H2 = (N_cells * I) / (2 * F) * 0.002016;
specific_energy = (P_stack / 1000) ./ (mass_H2 * 3600);
efficiency = 33.3 ./ specific_energy * 100;
energy_per_vol = specific_energy * 0.0899;

% Display at J = 1.0 A/cm²
idx_1A = find(j_range >= 1.0, 1, 'first');
fprintf('\n========================================\n');
fprintf('ENERGY ANALYSIS AT J = 1.0 A/cm²\n');
fprintf('========================================\n');
fprintf('  Specific energy: %.2f kWh/kg\n', specific_energy(idx_1A));
fprintf('  LHV Efficiency: %.2f %%\n', efficiency(idx_1A));
fprintf('  Stack Power: %.2f kW\n', P_stack(idx_1A)/1000);

% Create figure
fig = figure('Position', [100, 100, 1000, 700], 'Color', 'white');

% Subplot 1: Power
subplot(2, 2, 1);
plot(j_range, P_stack/1000, 'b-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Stack Power [kW]', 'FontSize', 10);
title('(1) Electrical Power Consumption', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);

% Subplot 2: Specific Energy
subplot(2, 2, 2);
plot(j_range, specific_energy, 'r-', 'LineWidth', 2.5);
hold on;
yline(33.3, 'k--', 'LHV (33.3 kWh/kg)', 'LineWidth', 1.5);
yline(39.4, 'k:', 'Thermodynamic Minimum (39.4 kWh/kg)', 'LineWidth', 1.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Specific Energy [kWh/kg H₂]', 'FontSize', 10);
title('(2) Specific Energy Consumption', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 8);
grid on; box on;
xlim([0, 2.0]);

% Subplot 3: Efficiency
subplot(2, 2, 3);
plot(j_range, efficiency, 'g-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('LHV Efficiency [%]', 'FontSize', 10);
title('(3) LHV-Based Energy Efficiency', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);

% Subplot 4: Energy per volume
subplot(2, 2, 4);
plot(j_range, energy_per_vol, 'm-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Energy per Volume [kWh/Nm³ H₂]', 'FontSize', 10);
title('(4) Energy Use per Unit Volume', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);

sgtitle('Energy Analysis of PEM Electrolyzer', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig, 'Energy_Analysis_of_PEM_Electrolyzer.jpg');
fprintf('\n✅ Figure saved: Energy_Analysis_of_PEM_Electrolyzer.jpg\n');
