%% FARADAY PRODUCTION ANALYSIS
% PEM Electrolyzer Faraday-Based Production Plots
% Uses validated 4-parameter model with literature-derived bounds
% k_conc = 4e-6 (literature-bounded)

clear; clc; close all;

%% ========== LOAD EXPERIMENTAL DATA ==========
% Update this path to where your PEMData.mat is located
loads = load('C:\Users\Admin\Desktop\Electrolyzer-modelling-toolbox-master\electrolyzerModellingToolbox\DemoData\PEMData.mat');
UIData = loads.UIData;
area = 17.64; % cm² (only used for validation, not needed for production curves)

%% ========== MODEL PARAMETERS (4-PARAMETER, LITERATURE BOUNDS) ==========
% Calibrated using fmincon with k_conc ∈ [0, 0.02]
% Combined RMSE = 0.018986 V, R² = 0.986711
E_rev = 1.42;               % Reversible potential at 75°C [V]
alpha = 0.448938;           % Transfer coefficient [-]
j0 = 5.428292e-02;          % Exchange current density [A/cm²]
R_ohm = 0.169475;           % Area-specific resistance [Ω·cm²]
k_conc = 4e-6;              % Concentration coefficient [V/(A/cm²)²]
n_e = 2;                    % Number of electrons

R = 8.314;                  % Gas constant [J/(mol·K)]
F = 96485;                  % Faraday constant [C/mol]
T = 348.15;                 % 75°C [K]

%% ========== MODEL FUNCTION ==========
model = @(j) E_rev + ...
            (R * T / (alpha * n_e * F)) * asinh(j / (2 * j0)) + ...
            R_ohm * j + k_conc * j.^2;

%% ========== POLARIZATION CURVE RANGE ==========
j_range = linspace(0.005, 2.2, 500);
V_range = model(j_range);

%% ========== FARADAY CALCULATIONS ==========
N_cells = 80;        % Representative stack size
A_cell = 424;        % cm² (Pareto‑optimal value from balanced solution)
I_f = j_range * A_cell;                     % Stack current [A]

n_H2 = (N_cells * I_f) / (2 * F);          % H₂ production [mol/s]
mass_H2_gph = n_H2 * 0.002016 * 3600 * 1000; % H₂ [g/h]

n_O2 = n_H2 / 2;                            % O₂ production [mol/s]
mass_O2_gph = n_O2 * 0.032 * 3600 * 1000;   % O₂ [g/h]
vol_O2_Lpm = n_O2 * 22.414 * 60;            % O₂ flow rate [L/min]

% Cumulative production over 24 hours
cum_H2_kg = n_H2 * 0.002016 * 3600 * 24;    % H₂ [kg/day]
cum_O2_kg = n_O2 * 0.032 * 3600 * 24;       % O₂ [kg/day]

%% ========== FIGURE 4: FARADAY PRODUCTION (4 subplots) ==========
fig4 = figure('Position', [100, 100, 1200, 700], 'Color', 'white');

% Subplot 1: H₂ and O₂ production rates (g/h)
subplot(2, 2, 1);
hold on; box on; grid on;
plot(j_range, mass_H2_gph, 'b-', 'LineWidth', 2.5, 'DisplayName', 'H₂ (g/h)');
plot(j_range, mass_O2_gph, 'r-', 'LineWidth', 2.5, 'DisplayName', 'O₂ (g/h)');
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('Production Rate [g/h]', 'FontSize', 10);
title('(a) H₂ and O₂ Production Rates', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 8);
xlim([0, 2.0]);
ylim([0, max(mass_H2_gph) * 1.1]);

% Subplot 2: O₂ flow rate (L/min) with 60 L/min target marked
subplot(2, 2, 2);
plot(j_range, vol_O2_Lpm, 'g-', 'LineWidth', 2.5);
hold on;
% Find and mark the 60 L/min target
idx_60 = find(vol_O2_Lpm >= 60, 1, 'first');
if ~isempty(idx_60)
    plot(j_range(idx_60), vol_O2_Lpm(idx_60), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(j_range(idx_60)+0.03, vol_O2_Lpm(idx_60)+3, ...
        sprintf('60 L/min at J = %.3f A/cm²', j_range(idx_60)), ...
        'FontSize', 9, 'FontWeight', 'bold');
end
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('O₂ Flow Rate [L/min]', 'FontSize', 10);
title('(b) Oxygen Flow Rate', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);
ylim([0, max(vol_O2_Lpm) * 1.1]);

% Subplot 3: Cumulative H₂ production (kg/day)
subplot(2, 2, 3);
plot(j_range, cum_H2_kg, 'b-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('H₂ Production [kg/day]', 'FontSize', 10);
title('(c) Cumulative H₂ Production (24 hours)', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);
ylim([0, max(cum_H2_kg) * 1.1]);

% Subplot 4: Cumulative O₂ production (kg/day)
subplot(2, 2, 4);
plot(j_range, cum_O2_kg, 'r-', 'LineWidth', 2.5);
xlabel('Current Density [A/cm²]', 'FontSize', 10);
ylabel('O₂ Production [kg/day]', 'FontSize', 10);
title('(d) Cumulative O₂ Production (24 hours)', 'FontSize', 11, 'FontWeight', 'bold');
grid on; box on;
xlim([0, 2.0]);
ylim([0, max(cum_O2_kg) * 1.1]);

% Main title
sgtitle('Faraday-Based Production Analysis (4-Parameter Model, k_{conc} = 4×10^{-6})', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(fig4, 'PEM_Faraday_Production.png');
fprintf('✅ Figure saved as PEM_Faraday_Production.png\n');

% Also save as PDF
print('PEM_Faraday_Production', '-dpdf', '-bestfit');
fprintf('✅ PDF saved as PEM_Faraday_Production.pdf\n');
