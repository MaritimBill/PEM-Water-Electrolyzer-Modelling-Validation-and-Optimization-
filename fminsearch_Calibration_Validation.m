%% FIX k_conc = 0 - PHYSICALLY DEFENSIBLE CALIBRATION
% This script calibrates the PEM electrolyzer model with k_conc = 0
% and generates a publication-ready validation figure.
% Colors: Black (Dataset 1), Red (Dataset 2), Dark Blue (Dataset 3)

clear; clc; close all;

fprintf('========================================\n');
fprintf('  PEM CALIBRATION: k_conc = 0\n');
fprintf('  Physically Defensible Model\n');
fprintf('========================================\n\n');

% Load data
loads = load('C:\Users\Admin\Desktop\Electrolyzer-modelling-toolbox-master\electrolyzerModellingToolbox\DemoData\PEMData.mat');
UIData = loads.UIData;
area = 17.64; % cm²

% Display data info
fprintf('📁 Data loaded:\n');
fprintf('  Electrode area: %.2f cm²\n', area);
fprintf('  Dataset 1: %d points (calibration)\n', length(UIData(1).I));
fprintf('  Dataset 2: %d points (validation)\n', length(UIData(2).I));
fprintf('  Dataset 3: %d points (validation)\n', length(UIData(3).I));
fprintf('  Total: %d points\n\n', length(UIData(1).I) + length(UIData(2).I) + length(UIData(3).I));

% Calibration data (Dataset 1)
U_exp = UIData(1).U;
I_exp = UIData(1).I;
j_exp = I_exp / area;

% Physical constants
T = 348.15; % 75°C
E_rev = 1.42;
R = 8.314;
F = 96485;
n = 2;

% Model with k_conc FIXED to 0
model = @(p, j) E_rev + ...
                (R * T / (p(1) * n * F)) * asinh(j / (2 * p(2))) + ...
                p(3) * j;  % NO k_conc term

% Objective: minimize RMSE
objfun = @(p) sqrt(mean((model(p, j_exp) - U_exp).^2));

% Initial guess: [alpha, j0, R_ohm]
p0 = [0.3, 1e-5, 0.18];

% Run fminsearch (unconstrained for the remaining 3 params)
fprintf('🔄 Running fminsearch optimization...\n\n');
options = optimset('Display', 'iter', 'MaxIter', 10000, 'TolX', 1e-12, 'TolFun', 1e-12);
[p_opt, fval, exitflag, output] = fminsearch(objfun, p0, options);

% Extract calibrated parameters
alpha = p_opt(1);
j0 = p_opt(2);
R_ohm = p_opt(3);
k_conc = 0;  % FIXED to zero

fprintf('\n========================================\n');
fprintf('  CALIBRATED PARAMETERS (k_conc = 0)\n');
fprintf('========================================\n');
fprintf('  Exit flag: %d\n', exitflag);
fprintf('  Iterations: %d\n', output.iterations);
fprintf('  Function evaluations: %d\n', output.funcCount);
fprintf('\n');
fprintf('  ┌─────────────────────────────────────────────┐\n');
fprintf('  │  CALIBRATED BUTLER-VOLMER PARAMETERS        │\n');
fprintf('  ├─────────────────────────────────────────────┤\n');
fprintf('  │  alpha  = %.6f  (transfer coefficient)      │\n', alpha);
fprintf('  │  j0     = %.6e A/cm²  (exchange current)    │\n', j0);
fprintf('  │  R_ohm  = %.6f Ω·cm²  (area resistance)     │\n', R_ohm);
fprintf('  │  k_conc = %.6f V/(A/cm²)²  (FIXED to zero)  │\n', k_conc);
fprintf('  │  RMSE   = %.6f V  (calibration dataset)     │\n', sqrt(fval / length(j_exp)));
fprintf('  └─────────────────────────────────────────────┘\n');

% Validate on all datasets
fprintf('\n========================================\n');
fprintf('  VALIDATION ON ALL DATASETS\n');
fprintf('========================================\n');

% Store results for plotting
all_U_exp = [];
all_U_pred = [];
all_rmse = zeros(3, 1);
all_mape = zeros(3, 1);
all_r2 = zeros(3, 1);
all_mean_error = zeros(3, 1);

for i = 1:3
    U_exp_i = UIData(i).U;
    I_exp_i = UIData(i).I;
    j_exp_i = I_exp_i / area;
    U_pred_i = model(p_opt, j_exp_i);
    
    % Metrics
    rmse_i = sqrt(mean((U_pred_i - U_exp_i).^2));
    mape_i = 100 * mean(abs((U_pred_i - U_exp_i) ./ U_exp_i));
    mean_error_i = mean(U_pred_i - U_exp_i);
    ss_res = sum((U_exp_i - U_pred_i).^2);
    ss_tot = sum((U_exp_i - mean(U_exp_i)).^2);
    r2_i = 1 - ss_res / ss_tot;
    
    all_rmse(i) = rmse_i;
    all_mape(i) = mape_i;
    all_r2(i) = r2_i;
    all_mean_error(i) = mean_error_i;
    
    fprintf('\nDataset %d:\n', i);
    fprintf('  RMSE = %.6f V\n', rmse_i);
    fprintf('  MAPE = %.4f %%\n', mape_i);
    fprintf('  R²   = %.6f\n', r2_i);
    fprintf('  Mean Error = %.6f V\n', mean_error_i);
    
    all_U_exp = [all_U_exp; U_exp_i];
    all_U_pred = [all_U_pred; U_pred_i];
end

% Combined metrics
combined_rmse = sqrt(mean((all_U_pred - all_U_exp).^2));
combined_mape = 100 * mean(abs((all_U_pred - all_U_exp) ./ all_U_exp));
combined_mean_error = mean(all_U_pred - all_U_exp);
ss_res = sum((all_U_exp - all_U_pred).^2);
ss_tot = sum((all_U_exp - mean(all_U_exp)).^2);
combined_r2 = 1 - ss_res / ss_tot;

fprintf('\n========================================\n');
fprintf('  COMBINED METRICS (All %d points)\n', length(all_U_exp));
fprintf('========================================\n');
fprintf('  RMSE = %.6f V\n', combined_rmse);
fprintf('  MAPE = %.4f %%\n', combined_mape);
fprintf('  R²   = %.6f\n', combined_r2);
fprintf('  Mean Error = %.6f V\n', combined_mean_error);

% Save parameters
params.alpha = alpha;
params.j0 = j0;
params.R_ohm = R_ohm;
params.k_conc = 0;
params.E_rev = E_rev;
params.area = area;
save('PEM_fit_params_k0.mat', 'params');
fprintf('\n✅ Parameters saved to: PEM_fit_params_k0.mat\n');

% Also save in PlatEMO format
p_opt = [alpha, j0, R_ohm, 0];
save('PEM_fit_params.mat', 'params', 'p_opt');
fprintf('✅ Parameters saved to: PEM_fit_params.mat (for PlatEMO)\n');

%% ========== PUBLICATION-READY PLOT ==========
% Colors: Black (Dataset 1), Red (Dataset 2), Dark Blue (Dataset 3)
colors = {
    [0, 0, 0],        % Black for Dataset 1
    [0.8, 0, 0],      % Red for Dataset 2
    [0, 0, 0.6]       % Dark Blue for Dataset 3
};
markers = {'o', 's', '^'};

figure('Color', 'white', 'Position', [100, 100, 1100, 500]);

%% Left Panel: Experimental vs Model Predictions
subplot(1, 2, 1);
hold on; box on; grid on;

for i = 1:3
    U_exp_i = UIData(i).U;
    I_exp_i = UIData(i).I;
    j_exp_i = I_exp_i / area;
    U_pred_i = model(p_opt, j_exp_i);
    
    % Experimental data (markers)
    plot(j_exp_i, U_exp_i, markers{i}, ...
        'Color', colors{i}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors{i}, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Dataset %d - Experimental', i));
    
    % Model predictions (lines)
    plot(j_exp_i, U_pred_i, '-', ...
        'Color', colors{i}, ...
        'LineWidth', 2.5, ...
        'DisplayName', sprintf('Dataset %d - Model', i));
end

xlabel('Current Density [A/cm²]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Cell Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold');
title('Validation: Experimental vs Model Predictions', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
xlim([0, 2.2]);
ylim([1.3, 2.1]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

%% Right Panel: Validation Metrics (RMSE and MAPE)
subplot(1, 2, 2);
hold on; box on; grid on;

% RMSE bars (left y-axis)
yyaxis left;
for i = 1:3
    bar(i, all_rmse(i), 0.4, 'FaceColor', colors{i}, 'EdgeColor', 'k', 'LineWidth', 1.5);
end
ylabel('RMSE [V]', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, 0.05]);
set(gca, 'YColor', 'k', 'FontSize', 12);

% MAPE markers (right y-axis)
yyaxis right;
for i = 1:3
    plot(i, all_mape(i), 'o', 'Color', colors{i}, 'MarkerSize', 10, 'LineWidth', 2);
end
plot(1:3, all_mape, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
ylabel('MAPE [%]', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, 2.5]);
set(gca, 'YColor', 'k', 'FontSize', 12);

xlabel('Dataset', 'FontSize', 14, 'FontWeight', 'bold');
title('Validation Metrics', 'FontSize', 14, 'FontWeight', 'bold');
xticks(1:3);
xticklabels({'Dataset 1', 'Dataset 2', 'Dataset 3'});

% Add combined metrics text box
metrics_text = sprintf(['Combined Metrics (59 points):\n' ...
                        'RMSE = %.6f V\n' ...
                        'MAPE = %.4f %%\n' ...
                        'R²   = %.6f'], ...
                        combined_rmse, combined_mape, combined_r2);

annotation('textbox', [0.54, 0.02, 0.15, 0.20], ...
    'String', metrics_text, ...
    'FontSize', 10, ...
    'BackgroundColor', [0.95, 0.95, 0.95], ...
    'EdgeColor', [0.5, 0.5, 0.5], ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center');

%% Add legend for right panel (RMSE and MAPE)
% Create invisible axes for legend
ax_legend = axes('Position', [0.88, 0.70, 0.10, 0.15], 'Visible', 'off');
hold on;
bar(1, 0.01, 0.4, 'FaceColor', 'k', 'EdgeColor', 'k');
plot(1, 0.01, 'o', 'Color', 'k', 'MarkerSize', 10);
text(0.5, 0.02, 'RMSE (bars)', 'FontSize', 9, 'HorizontalAlignment', 'right');
text(0.5, 0.01, 'MAPE (markers)', 'FontSize', 9, 'HorizontalAlignment', 'right');
axis off;

%% Main Title
sgtitle(sprintf('PEM Electrolyzer Validation (k_{conc} = 0): RMSE = %.4f V, R² = %.6f', ...
    combined_rmse, combined_r2), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Save Figure
saveas(gcf, 'PEM_Current_Validation_k0.png');
fprintf('✅ Figure saved to: PEM_Current_Validation_k0.png\n');

%% Save as high-resolution PDF
set(gcf, 'PaperPositionMode', 'auto');
print('PEM_Current_Validation_k0', '-dpdf', '-bestfit');
fprintf('✅ PDF saved to: PEM_Current_Validation_k0.pdf\n');

%% ========== FINAL SUMMARY ==========
fprintf('\n========================================\n');
fprintf('  CALIBRATION COMPLETE\n');
fprintf('========================================\n');
fprintf('\n📊 Summary of Results (k_conc = 0):\n');
fprintf('  ┌─────────────────────────────────────────────────────┐\n');
fprintf('  │ Model: Butler-Volmer with k_conc = 0              │\n');
fprintf('  │ α     = %.6f                                      │\n', alpha);
fprintf('  │ j₀    = %.6e A/cm²                               │\n', j0);
fprintf('  │ R_ohm = %.6f Ω·cm²                               │\n', R_ohm);
fprintf('  │ k_conc= %.6f V/(A/cm²)² (FIXED)                 │\n', k_conc);
fprintf('  │ RMSE  = %.6f V (combined, %d points)             │\n', combined_rmse, length(all_U_exp));
fprintf('  │ R²    = %.6f                                     │\n', combined_r2);
fprintf('  └─────────────────────────────────────────────────────┘\n');

fprintf('\n✅ All done!\n');