%% ========== COMPLETE SCRIPT: VALIDATION + CSV EXPORT ==========
% Load parameters (if not already in workspace)
if ~exist('p_opt', 'var')
    if exist('PEM_fit_params.mat', 'file')
        load('PEM_fit_params.mat', 'p_opt');
    else
        error('Run calibration first to create PEM_fit_params.mat');
    end
end

% Load data - adjust path below to match your system!
loads = load('C:\Users\Admin\Desktop\Electrolyzer-modelling-toolbox-master\electrolyzerModellingToolbox\DemoData\PEMData.mat');
UIData = loads.UIData;
area = 17.64;

% Model definition (k_conc=0)
T = 348.15; E_rev = 1.42; R = 8.314; F = 96485; n = 2;
model = @(p, j) E_rev + (R*T/(p(1)*n*F))*asinh(j/(2*p(2))) + p(3)*j + p(4)*j.^2;

% Compute predictions and metrics
all_U_exp = []; all_U_pred = [];
all_rmse = zeros(3,1); all_mape = zeros(3,1);
for i = 1:3
    U_exp_i = UIData(i).U;
    I_exp_i = UIData(i).I;
    j_exp_i = I_exp_i / area;
    U_pred_i = model(p_opt, j_exp_i);
    all_U_exp = [all_U_exp; U_exp_i];
    all_U_pred = [all_U_pred; U_pred_i];
    all_rmse(i) = sqrt(mean((U_pred_i - U_exp_i).^2));
    all_mape(i) = 100 * mean(abs((U_pred_i - U_exp_i) ./ U_exp_i));
end
combined_rmse = sqrt(mean((all_U_pred - all_U_exp).^2));
combined_mape = 100 * mean(abs((all_U_pred - all_U_exp) ./ all_U_exp));
ss_res = sum((all_U_exp - all_U_pred).^2);
ss_tot = sum((all_U_exp - mean(all_U_exp)).^2);
combined_r2 = 1 - ss_res / ss_tot;

% Colors and marker style
colors = {[0,0,0], [0.8,0,0], [0,0,0.6]};  % Black, Red, Dark Blue
marker = 'o';   % circles only

% Create figure – new size: 1355 x 690 pixels
figure('Color', 'white', 'Position', [100, 100, 1355, 690]);

% Define custom positions: [left, bottom, width, height] in normalized units
% First subplot (I‑V) – wider: 38% of figure width
pos1 = [0.08, 0.15, 0.38, 0.75];
% Second subplot (RMSE) – 22% width
pos2 = [0.50, 0.15, 0.22, 0.75];
% Third subplot (MAPE) – 22% width
pos3 = [0.76, 0.15, 0.22, 0.75];

%% ----- Subplot 1: I‑V curves (wider) -----
ax1 = axes('Position', pos1);
hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
for i = 1:3
    U_exp_i = UIData(i).U;
    I_exp_i = UIData(i).I;
    j_exp_i = I_exp_i / area;
    U_pred_i = model(p_opt, j_exp_i);
    % Model line (hidden from legend)
    plot(ax1, j_exp_i, U_pred_i, '-', 'Color', colors{i}, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % Experimental markers (circles)
    h(i) = plot(ax1, j_exp_i, U_exp_i, marker, ...
        'Color', colors{i}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors{i}, ...
        'MarkerEdgeColor', colors{i}, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Dataset %d', i));
end
xlabel(ax1, 'Current Density [A/cm²]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax1, 'Cell Voltage [V]', 'FontSize', 14, 'FontWeight', 'bold');
title(ax1, '(a) Experimental vs Model', 'FontSize', 14, 'FontWeight', 'bold');
xlim(ax1, [0, 2.2]); ylim(ax1, [1.3, 2.1]);
legend(ax1, h, 'Location', 'southeast', 'FontSize', 12);
set(ax1, 'FontSize', 12, 'LineWidth', 1.5);

%% ----- Subplot 2: RMSE (bars) -----
ax2 = axes('Position', pos2);
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
b = bar(ax2, 1:3, all_rmse, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.5);
for i = 1:3
    b.CData(i,:) = colors{i};
    text(ax2, i, all_rmse(i) + 0.001, sprintf('%.4f', all_rmse(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
end
xlabel(ax2, 'Dataset', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax2, 'RMSE [V]', 'FontSize', 14, 'FontWeight', 'bold');
title(ax2, '(b) Root Mean Square Error', 'FontSize', 14, 'FontWeight', 'bold');
xticks(ax2, 1:3); xticklabels(ax2, {'1', '2', '3'});
ylim(ax2, [0, max(all_rmse)*1.3]);
set(ax2, 'FontSize', 12, 'LineWidth', 1.5);

%% ----- Subplot 3: MAPE (bars) -----
ax3 = axes('Position', pos3);
hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on');
b2 = bar(ax3, 1:3, all_mape, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.5);
for i = 1:3
    b2.CData(i,:) = colors{i};
    text(ax3, i, all_mape(i) + 0.02, sprintf('%.2f', all_mape(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
end
xlabel(ax3, 'Dataset', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax3, 'MAPE [%]', 'FontSize', 14, 'FontWeight', 'bold');
title(ax3, '(c) Mean Absolute Percentage Error', 'FontSize', 14, 'FontWeight', 'bold');
xticks(ax3, 1:3); xticklabels(ax3, {'1', '2', '3'});
ylim(ax3, [0, max(all_mape)*1.3]);
set(ax3, 'FontSize', 12, 'LineWidth', 1.5);

%% ----- Save figure (PNG and PDF) -----
saveas(gcf, 'PEM_Validation_3panels_wider.png');
fprintf('✅ Figure saved as PEM_Validation_3panels_wider.png\n');
print('PEM_Validation_3panels_wider', '-dpdf', '-bestfit');
fprintf('✅ PDF saved as PEM_Validation_3panels_wider.pdf\n');

%% ====== EXPORT DATA TO CSV ======
% Open the CSV file for writing
fid = fopen('PEM_Validation_Data.csv', 'w');
% Write the header row
fprintf(fid, 'Dataset,CurrentDensity_A_cm2,Voltage_Exp_V,Voltage_Pred_V,Error_V\n');

for i = 1:3
    % Retrieve data for the current dataset again
    U_exp_i = UIData(i).U;
    I_exp_i = UIData(i).I;
    j_exp_i = I_exp_i / area;
    U_pred_i = model(p_opt, j_exp_i);
    
    % Write each data point for this dataset (using 16 decimal places for high precision)
    for k = 1:length(j_exp_i)
        fprintf(fid, '%d,%.16f,%.16f,%.16f,%.16f\n', ...
            i, ...
            j_exp_i(k), ...
            U_exp_i(k), ...
            U_pred_i(k), ...
            U_exp_i(k) - U_pred_i(k));  % Calculates Error_V (Exp - Pred)
    end
end

% Close the file
fclose(fid);
fprintf('✅ CSV data saved as PEM_Validation_Data.csv\n');
