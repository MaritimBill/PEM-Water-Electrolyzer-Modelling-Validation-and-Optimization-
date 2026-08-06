%% ========================================================================
%% PARAMETER EFFECTS ON PEM POLARIZATION CURVES - FIGURE 5 (Paper)
%% Using Butler-Volmer model exactly as in paper
%% Fixed parameters: α = 0.4, j0 = 1.5×10⁻⁶ A/cm², jL = 1.5 A/cm², r = 1 Ω·cm²
%% ========================================================================
clear; clc; close all;

%% ============================================================
%% PARAMETERS (from paper - Figure 5 caption)
%% ============================================================
alpha_fixed = 0.4;
j0_fixed = 1.5e-6;
jL_fixed = 1.5;
r_fixed = 1.0;

% Temperature
T = 75 + 273.15;  % 75°C in Kelvin

% Constants
F = 96485;
R = 8.314;
n_e = 2;

% Reversible potential at 75°C
E_rev = 1.42;

% Current density range
j = linspace(0.001, 2.0, 300)';

%% Create figure
figure('Position', [100, 100, 1200, 900], 'Color', 'white');

% Colors for plots
colors = lines(3);

% ============================================================
% (a) Electron transfer coefficient (α)
% ============================================================
subplot(2, 2, 1);
hold on; grid on; box on;

alpha_vals = [0.25, 0.4, 0.75];
j0 = j0_fixed;
r = r_fixed;
jL = jL_fixed;

for i = 1:length(alpha_vals)
    V = calcVoltageVec(j, alpha_vals(i), j0, r, jL, E_rev, T);
    V(V > 5) = NaN;
    plot(j, V, 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('\\alpha = %.2f', alpha_vals(i)));
end

xlabel('Current Density, j (A/cm²)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cell Voltage, U (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('(a) Electron Transfer Coefficient (\alpha)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([0, 1.8]); ylim([0, 4.5]);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'LineWidth', 1.2);

% Add magnified inset
axes('Position', [0.55, 0.68, 0.25, 0.18]);
hold on;
for i = 1:length(alpha_vals)
    V = calcVoltageVec(j, alpha_vals(i), j0, r, jL, E_rev, T);
    V(V > 5) = NaN;
    plot(j, V, 'Color', colors(i,:), 'LineWidth', 2);
end
xlim([0, 0.3]); ylim([1.4, 1.55]);
grid on;
title('Low j region', 'FontSize', 9);
set(gca, 'FontSize', 8);

% ============================================================
% (b) Exchange current density (j0)
% ============================================================
subplot(2, 2, 2);
hold on; grid on; box on;

j0_vals = [1.5e-8, 1.5e-6, 1.5e-4];
alpha = alpha_fixed;
r = r_fixed;
jL = jL_fixed;

for i = 1:length(j0_vals)
    V = calcVoltageVec(j, alpha, j0_vals(i), r, jL, E_rev, T);
    V(V > 5) = NaN;
    plot(j, V, 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('j_0 = %.1e A/cm²', j0_vals(i)));
end

xlabel('Current Density, j (A/cm²)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cell Voltage, U (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) Exchange Current Density (j_0)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([0, 1.8]); ylim([0, 4.5]);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'LineWidth', 1.2);

% Add magnified inset
axes('Position', [0.55, 0.20, 0.25, 0.18]);
hold on;
for i = 1:length(j0_vals)
    V = calcVoltageVec(j, alpha, j0_vals(i), r, jL, E_rev, T);
    V(V > 5) = NaN;
    plot(j, V, 'Color', colors(i,:), 'LineWidth', 2);
end
xlim([0, 0.3]); ylim([1.4, 1.55]);
grid on;
title('Low j region', 'FontSize', 9);
set(gca, 'FontSize', 8);

% ============================================================
% (c) Limiting current density (jL)
% ============================================================
subplot(2, 2, 3);
hold on; grid on; box on;

jL_vals = [1.0, 1.5, 2.0];
alpha = alpha_fixed;
j0 = j0_fixed;
r = r_fixed;

for i = 1:length(jL_vals)
    V = calcVoltageVec(j, alpha, j0, r, jL_vals(i), E_rev, T);
    V(V > 5) = NaN;
    plot(j, V, 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('j_L = %.1f A/cm²', jL_vals(i)));
end

xlabel('Current Density, j (A/cm²)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cell Voltage, U (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('(c) Limiting Current Density (j_L)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([0, 1.8]); ylim([0, 4.5]);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'LineWidth', 1.2);

% ============================================================
% (d) Area-specific resistance (r)
% ============================================================
subplot(2, 2, 4);
hold on; grid on; box on;

r_vals = [0.5, 1.0, 1.5];
alpha = alpha_fixed;
j0 = j0_fixed;
jL = jL_fixed;

for i = 1:length(r_vals)
    V = calcVoltageVec(j, alpha, j0, r_vals(i), jL, E_rev, T);
    V(V > 5) = NaN;
    plot(j, V, 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('r = %.1f \\Omega cm²', r_vals(i)));
end

xlabel('Current Density, j (A/cm²)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cell Voltage, U (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('(d) Area-Specific Resistance (r)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([0, 1.8]); ylim([0, 4.5]);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'LineWidth', 1.2);

%% Main title
sgtitle('Effect of Different Fitting Parameters on the Polarization Curve', 'FontSize', 14, 'FontWeight', 'bold');

%% Save figure
saveas(gcf, 'PEM_Parameter_Effects_v4.png');
fprintf('Saved: PEM_Parameter_Effects_v4.png\n');
fprintf('Using Butler-Volmer model as in paper\n');
fprintf('Fixed params: alpha=%.2f, j0=%.1e, jL=%.1f, r=%.1f\n', alpha_fixed, j0_fixed, jL_fixed, r_fixed);

%% ============================================================
%% PEM POLARIZATION MODEL - Butler-Volmer (vectorized)
%% ============================================================
function V = calcVoltageVec(j, alpha, j0, r, jL, E_rev, T)
    F = 96485;
    R = 8.314;
    n_e = 2;
    
    % Activation overpotential - Butler-Volmer (vectorized)
    f = R * T / (alpha * n_e * F);
    eta_act = f * asinh(j / (2 * j0));
    eta_act(j <= 0) = 0;
    
    % Ohmic overpotential
    eta_ohm = r * j;
    
    % Concentration overpotential (vectorized)
    eta_conc = zeros(size(j));
    valid_idx = (j < jL) & (jL > 0);
    eta_conc(valid_idx) = (R * T / (n_e * F)) * log(jL ./ (jL - j(valid_idx)));
    eta_conc(~valid_idx) = 10;  % Large penalty for exceeding limit
    
    % Total voltage
    V = E_rev + eta_act + eta_ohm + eta_conc;
end
