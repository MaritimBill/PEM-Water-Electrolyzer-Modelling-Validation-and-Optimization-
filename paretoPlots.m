%% Improved Pareto Front Visualization - LAST RUN ONLY (With Solution Table)
clear; clc; close all;

% Add PlatEMO to path
addpath(genpath('C:\Users\Admin\Desktop\PlatEMO-4.7'));

%% Problem Definition
problem = @PEM_Problem_Scaled;
popSize = 200;
maxGen = 100;
maxFE = popSize * maxGen;
numRuns = 30;

fprintf('========================================\n');
fprintf('Pareto Visualization (Last Run Only)\n');
fprintf('========================================\n');

%% Run FeaNSGA_II
final_objs = [];
final_decs = [];
final_cons = [];

for run = 1:numRuns
    fprintf('\n🔹 Run %d/%d...\n', run, numRuns);
    
    Problem = problem();
    Problem.N = popSize;
    Problem.maxFE = maxFE;
    
    Algorithm = FeaNSGA_II();
    tic;
    Algorithm.Solve(Problem);
    elapsed = toc;
    
    if ~isempty(Algorithm.result)
        finalPop = Algorithm.result{end,2};
        
        n = length(finalPop);
        objs = zeros(n, 3);
        decs = zeros(n, 5);
        cons = zeros(n, 10);
        
        for i = 1:n
            objs_i = finalPop(i).objs;
            if size(objs_i,1) > 1, objs_i = objs_i'; end
            objs(i,:) = objs_i;
            
            decs_i = finalPop(i).decs;
            if size(decs_i,1) > 1, decs_i = decs_i'; end
            decs(i,:) = decs_i;
            
            cons_i = finalPop(i).cons;
            if size(cons_i,1) > 1, cons_i = cons_i'; end
            cons(i,:) = cons_i;
        end
        
        feasible = all(cons <= 0, 2);
        
        if any(feasible)
            objs_feas = objs(feasible,:);
            decs_feas = decs(feasible,:);
            
            ND = NDSort(objs_feas, 1) == 1;
            pareto_size = sum(ND);
            
            fprintf('   ✅ Feasible: %d, Pareto: %d, Time: %.2fs\n', ...
                sum(feasible), pareto_size, elapsed);
            
            final_objs = objs_feas;
            final_decs = decs_feas;
            final_cons = cons(feasible,:);
        end
    end
end

%% Get Pareto front of final run
if isempty(final_objs), error('No feasible solutions!'); end

ND_final = NDSort(final_objs, 1) == 1;
pareto_final = final_objs(ND_final, :);
pareto_decs_final = final_decs(ND_final, :);

fprintf('\n📊 FINAL RUN RESULTS:\n');
fprintf('   Total feasible: %d\n', size(final_objs,1));
fprintf('   Pareto solutions: %d (%.1f%% of population)\n', ...
    size(pareto_final,1), size(pareto_final,1)/popSize*100);

%% ========== CALCULATE HYPERVOLUME CORRECTLY ==========
% Normalize objectives to [0,1] range
minF = min(pareto_final, [], 1);
maxF = max(pareto_final, [], 1);
pareto_norm = (pareto_final - minF) ./ (maxF - minF + 1e-12);

% Reference point (worse than all normalized solutions)
ref_norm = [1.1, 1.1, 1.1];

% Manual HV calculation (since PlatEMO's Hypervolume might be picky)
try
    % Method 1: Use PlatEMO's Hypervolume
    hv_value = Hypervolume(pareto_norm, ref_norm);
    if isnan(hv_value) || hv_value <= 0
        error('PlatEMO HV failed');
    end
catch
    % Method 2: Manual Monte Carlo HV calculation
    N_samples = 10000;
    samples = rand(N_samples, 3) .* ref_norm;
    dominated = false(N_samples, 1);
    
    for i = 1:size(pareto_norm, 1)
        is_dominated = all(samples >= pareto_norm(i, :), 2);
        dominated = dominated | is_dominated;
    end
    
    hv_value = sum(dominated) / N_samples * prod(ref_norm);
end

fprintf('\n📊 HYPERVOLUME (normalized, ref=[1.1,1.1,1.1]): %.4f\n', hv_value);

%% ========== DISPLAY SOLUTION TABLE ==========
fprintf('\n%s\n', repmat('=', 1, 100));
fprintf('PARETO SOLUTIONS (%d total)\n', size(pareto_final,1));
fprintf('%s\n', repmat('=', 1, 100));
fprintf('%-6s %12s %12s %12s | %12s %12s %12s %12s %12s\n', ...
    'Index', 'Energy', 'Cost', 'Stress', 'A_cell', 'N_cells', 'J', 'T_C', 'P(bar)');
fprintf('%s\n', repmat('-', 1, 110));

nShow = min(50, size(pareto_final,1));
step = max(1, floor(size(pareto_final,1)/nShow));

for i = 1:step:size(pareto_final,1)
    fprintf('%-6d %12.4f %12.4f %12.4f | %12.2f %12.0f %12.0f %12.1f %12.2f\n', ...
        i, ...
        pareto_final(i,1), ...
        pareto_final(i,2), ...
        pareto_final(i,3), ...
        pareto_decs_final(i,1), ...
        pareto_decs_final(i,2), ...
        pareto_decs_final(i,3), ...
        pareto_decs_final(i,4), ...
        pareto_decs_final(i,5)/1e5);
end

if size(pareto_final,1) > 50
    fprintf('\n... and %d more solutions\n', size(pareto_final,1) - 50);
end
fprintf('%s\n', repmat('=', 1, 100));

%% ========== PUBLICATION QUALITY FIGURES ==========

set(0, 'DefaultAxesFontSize', 11, 'DefaultTextFontSize', 11);

%% Figure 1: 3D Pareto Front (Final Run Only)
fig1 = figure('Position', [50 50 1000 800], 'Color', 'white');

% Normalize for visualization
pareto_norm_vis = zeros(size(pareto_final));
for i = 1:3
    pareto_norm_vis(:,i) = (pareto_final(:,i) - min(pareto_final(:,i))) / ...
                       (max(pareto_final(:,i)) - min(pareto_final(:,i)) + 1e-10);
end

scatter3(pareto_norm_vis(:,1), pareto_norm_vis(:,2), pareto_norm_vis(:,3), 100, ...
    pareto_norm_vis(:,3), 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8);
colormap(jet);
colorbar;
xlabel('Energy (normalized)', 'FontWeight', 'bold');
ylabel('Cost (normalized)', 'FontWeight', 'bold');
zlabel('Stress (normalized)', 'FontWeight', 'bold');
title(sprintf('3D Pareto Front - %d Solutions (HV=%.4f)', size(pareto_final,1), hv_value), ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on; box on;
view(45, 30);
rotate3d on;

%% Figure 2: 2D Projections
fig2 = figure('Position', [100 100 1400 500], 'Color', 'white');

% Energy vs Cost
subplot(1,3,1);
scatter(pareto_final(:,1), pareto_final(:,2), 80, pareto_final(:,3), ...
    'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
colormap(jet); 
colorbar;
xlabel('Energy (kWh/kg)', 'FontWeight', 'bold');
ylabel('Cost (m²)', 'FontWeight', 'bold');
title(sprintf('Energy vs Cost (HV=%.4f)', hv_value), 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% Energy vs Stress
subplot(1,3,2);
scatter(pareto_final(:,1), pareto_final(:,3), 80, pareto_final(:,2), ...
    'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
colormap(jet); 
colorbar;
xlabel('Energy (kWh/kg)', 'FontWeight', 'bold');
ylabel('Stress Index', 'FontWeight', 'bold');
title('Energy vs Stress', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% Cost vs Stress
subplot(1,3,3);
scatter(pareto_final(:,2), pareto_final(:,3), 80, pareto_final(:,1), ...
    'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
colormap(jet); 
colorbar;
xlabel('Cost (m²)', 'FontWeight', 'bold');
ylabel('Stress Index', 'FontWeight', 'bold');
title('Cost vs Stress', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

sgtitle(sprintf('FeaNSGA_II: Pareto Front (%d solutions, HV=%.4f)', size(pareto_final,1), hv_value), ...
    'FontSize', 15, 'FontWeight', 'bold');

%% Figure 3: Decision Variable Distribution
fig3 = figure('Position', [150 150 1200 400], 'Color', 'white');

varNames = {'A_{cell} (cm²)', 'N_{cells}', 'J (A/m²)', 'T (°C)', 'P (bar)'};

for i = 1:5
    subplot(1,5,i);
    if i == 5
        values = pareto_decs_final(:,i) / 1e5;
        xlabel('P (bar)', 'FontWeight', 'bold');
    else
        values = pareto_decs_final(:,i);
        xlabel(varNames{i}, 'FontWeight', 'bold');
    end
    histogram(values, 15, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'k');
    if i == 1, ylabel('Count'); end
    title(sprintf('Range: %.2f-%.2f', min(values), max(values)), 'FontSize', 10);
    grid on; box on;
end
sgtitle('Pareto Decision Variable Distributions', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 4: Correlation Heatmap
fig4 = figure('Position', [200 200 600 500], 'Color', 'white');

corr_objs = corrcoef(pareto_final);
corr_decs = corrcoef(pareto_decs_final);

subplot(1,2,1);
imagesc(corr_objs);
colorbar;
set(gca, 'XTick', 1:3, 'YTick', 1:3, 'XTickLabel', {'Energy','Cost','Stress'}, ...
    'YTickLabel', {'Energy','Cost','Stress'});
title('Objective Correlations', 'FontWeight', 'bold');
colormap(jet);
axis square;

subplot(1,2,2);
imagesc(corr_decs);
colorbar;
set(gca, 'XTick', 1:5, 'YTick', 1:5, 'XTickLabel', {'Ac','Nc','J','T','P'}, ...
    'YTickLabel', {'Ac','Nc','J','T','P'});
title('Decision Variable Correlations', 'FontWeight', 'bold');
colormap(jet);
axis square;

%% ========== CONSOLE OUTPUT ==========
fprintf('\n📈 STATISTICAL SUMMARY (Final Run)\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('Objectives:\n');
fprintf('  Energy: %.4f - %.4f kWh/kg\n', min(pareto_final(:,1)), max(pareto_final(:,1)));
fprintf('  Cost:   %.4f - %.4f m²\n', min(pareto_final(:,2)), max(pareto_final(:,2)));
fprintf('  Stress: %.4f - %.4f\n', min(pareto_final(:,3)), max(pareto_final(:,3)));

fprintf('\nDecision Variables:\n');
fprintf('  A_cell:   %.2f - %.2f cm²\n', min(pareto_decs_final(:,1)), max(pareto_decs_final(:,1)));
fprintf('  N_cells:  %.0f - %.0f\n', min(pareto_decs_final(:,2)), max(pareto_decs_final(:,2)));
fprintf('  J:        %.0f - %.0f A/m²\n', min(pareto_decs_final(:,3)), max(pareto_decs_final(:,3)));
fprintf('  T_C:      %.1f - %.1f °C\n', min(pareto_decs_final(:,4)), max(pareto_decs_final(:,4)));
fprintf('  P:        %.2f - %.2f bar\n', min(pareto_decs_final(:,5))/1e5, max(pareto_decs_final(:,5))/1e5);

%% Best Solutions
fprintf('\n🏆 BEST SOLUTIONS (Final Run)\n');
[~, idx1] = min(pareto_final(:,1));
fprintf('  Min Energy:  Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(idx1,:));
[~, idx2] = min(pareto_final(:,2));
fprintf('  Min Cost:    Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(idx2,:));
[~, idx3] = min(pareto_final(:,3));
fprintf('  Min Stress:  Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(idx3,:));

%% Save figures
try
    exportgraphics(fig1, 'pareto_3D.png', 'Resolution', 300);
    exportgraphics(fig2, 'pareto_2D.png', 'Resolution', 300);
    exportgraphics(fig3, 'pareto_dec_vars.png', 'Resolution', 300);
    exportgraphics(fig4, 'pareto_correlations.png', 'Resolution', 300);
    fprintf('\n✅ Figures saved as PNG (300 DPI)\n');
catch
    saveas(fig1, 'pareto_3D.png');
    saveas(fig2, 'pareto_2D.png');
    saveas(fig3, 'pareto_dec_vars.png');
    saveas(fig4, 'pareto_correlations.png');
    fprintf('\n✅ Figures saved as PNG\n');
end

% Save data
save('pareto_analysis.mat', 'pareto_final', 'pareto_decs_final', 'final_objs', 'final_decs', 'hv_value');

fprintf('✅ Data saved: pareto_analysis.mat\n');
fprintf('✅ Hypervolume: %.4f\n', hv_value);
fprintf('\n%s\n', repmat('=', 1, 50));
fprintf('VISUALIZATION COMPLETE!\n');

%% ========== PRINT NUMERICAL VALUES FOR FIGURE CAPTIONS ==========
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('📊 EXACT VALUES FOR FIGURE CAPTIONS\n');
fprintf('%s\n', repmat('=', 1, 60));

% Get exact ranges
A_cell_min = min(pareto_decs_final(:,1));
A_cell_max = max(pareto_decs_final(:,1));
N_cells_min = min(pareto_decs_final(:,2));
N_cells_max = max(pareto_decs_final(:,2));
J_min = min(pareto_decs_final(:,3));
J_max = max(pareto_decs_final(:,3));
T_min = min(pareto_decs_final(:,4));
T_max = max(pareto_decs_final(:,4));
P_min = min(pareto_decs_final(:,5))/1e5;
P_max = max(pareto_decs_final(:,5))/1e5;

Energy_min = min(pareto_final(:,1));
Energy_max = max(pareto_final(:,1));
Cost_min = min(pareto_final(:,2));
Cost_max = max(pareto_final(:,2));
Stress_min = min(pareto_final(:,3));
Stress_max = max(pareto_final(:,3));

% Print Figure 3 values (Decision Variables)
fprintf('\n📊 FIGURE 3: Decision Variable Distributions\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('  A_cell:   %.2f - %.2f cm²\n', A_cell_min, A_cell_max);
fprintf('  N_cells:  %.0f - %.0f\n', N_cells_min, N_cells_max);
fprintf('  J:        %.0f - %.0f A/m²\n', J_min, J_max);
fprintf('  T_C:      %.1f - %.1f °C\n', T_min, T_max);
fprintf('  P:        %.2f - %.2f bar\n', P_min, P_max);

% Print Figure 1/2 values (Objectives)
fprintf('\n📊 FIGURE 1 & 2: Pareto Front Objectives\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('  Energy:   %.4f - %.4f kWh/kg\n', Energy_min, Energy_max);
fprintf('  Cost:     %.4f - %.4f m²\n', Cost_min, Cost_max);
fprintf('  Stress:   %.4f - %.4f\n', Stress_min, Stress_max);

% Print correlation matrices
corr_objs = corrcoef(pareto_final);
corr_decs = corrcoef(pareto_decs_final);

fprintf('\n📊 FIGURE 4: Correlation Matrices\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('Objective Correlations (Energy, Cost, Stress):\n');
fprintf('  Energy-Cost:   %.3f\n', corr_objs(1,2));
fprintf('  Energy-Stress: %.3f\n', corr_objs(1,3));
fprintf('  Cost-Stress:   %.3f\n', corr_objs(2,3));

fprintf('\nDecision Variable Correlations (A_cell, N_cells, J, T, P):\n');
fprintf('  A_cell-N_cells: %.3f\n', corr_decs(1,2));
fprintf('  A_cell-J:       %.3f\n', corr_decs(1,3));
fprintf('  A_cell-T:       %.3f\n', corr_decs(1,4));
fprintf('  A_cell-P:       %.3f\n', corr_decs(1,5));
fprintf('  N_cells-J:      %.3f\n', corr_decs(2,3));
fprintf('  N_cells-T:      %.3f\n', corr_decs(2,4));
fprintf('  N_cells-P:      %.3f\n', corr_decs(2,5));
fprintf('  J-T:            %.3f\n', corr_decs(3,4));
fprintf('  J-P:            %.3f\n', corr_decs(3,5));
fprintf('  T-P:            %.3f\n', corr_decs(4,5));

% Print sample solutions for Table 6
fprintf('\n📊 TABLE 6: Representative Pareto-Optimal Solutions\n');
fprintf('%s\n', repmat('-', 1, 50));

[~, minEnergyIdx] = min(pareto_final(:,1));
[~, minCostIdx] = min(pareto_final(:,2));
[~, minStressIdx] = min(pareto_final(:,3));

% Find a balanced solution (closest to ideal point)
ideal = [min(pareto_final(:,1)), min(pareto_final(:,2)), min(pareto_final(:,3))];
distances = sqrt(sum((pareto_final - ideal).^2, 2));
[~, balancedIdx] = min(distances);

fprintf('Min Energy:  Energy=%.2f, Cost=%.3f, Stress=%.2f | A_cell=%.0f, N_cells=%.0f, J=%.0f, T_C=%.1f, P=%.2f\n', ...
    pareto_final(minEnergyIdx,1), pareto_final(minEnergyIdx,2), pareto_final(minEnergyIdx,3), ...
    pareto_decs_final(minEnergyIdx,1), pareto_decs_final(minEnergyIdx,2), pareto_decs_final(minEnergyIdx,3), ...
    pareto_decs_final(minEnergyIdx,4), pareto_decs_final(minEnergyIdx,5)/1e5);

fprintf('Min Cost:    Energy=%.2f, Cost=%.3f, Stress=%.2f | A_cell=%.0f, N_cells=%.0f, J=%.0f, T_C=%.1f, P=%.2f\n', ...
    pareto_final(minCostIdx,1), pareto_final(minCostIdx,2), pareto_final(minCostIdx,3), ...
    pareto_decs_final(minCostIdx,1), pareto_decs_final(minCostIdx,2), pareto_decs_final(minCostIdx,3), ...
    pareto_decs_final(minCostIdx,4), pareto_decs_final(minCostIdx,5)/1e5);

fprintf('Min Stress:  Energy=%.2f, Cost=%.3f, Stress=%.2f | A_cell=%.0f, N_cells=%.0f, J=%.0f, T_C=%.1f, P=%.2f\n', ...
    pareto_final(minStressIdx,1), pareto_final(minStressIdx,2), pareto_final(minStressIdx,3), ...
    pareto_decs_final(minStressIdx,1), pareto_decs_final(minStressIdx,2), pareto_decs_final(minStressIdx,3), ...
    pareto_decs_final(minStressIdx,4), pareto_decs_final(minStressIdx,5)/1e5);

fprintf('Balanced:    Energy=%.2f, Cost=%.3f, Stress=%.2f | A_cell=%.0f, N_cells=%.0f, J=%.0f, T_C=%.1f, P=%.2f\n', ...
    pareto_final(balancedIdx,1), pareto_final(balancedIdx,2), pareto_final(balancedIdx,3), ...
    pareto_decs_final(balancedIdx,1), pareto_decs_final(balancedIdx,2), pareto_decs_final(balancedIdx,3), ...
    pareto_decs_final(balancedIdx,4), pareto_decs_final(balancedIdx,5)/1e5);

% Save to text file
fid = fopen('figure_values_for_captions.txt', 'w');
fprintf(fid, '========================================\n');
fprintf(fid, 'EXACT VALUES FOR FIGURE CAPTIONS\n');
fprintf(fid, '========================================\n\n');

fprintf(fid, 'FIGURE 3: Decision Variable Distributions\n');
fprintf(fid, '  A_cell:   %.2f - %.2f cm²\n', A_cell_min, A_cell_max);
fprintf(fid, '  N_cells:  %.0f - %.0f\n', N_cells_min, N_cells_max);
fprintf(fid, '  J:        %.0f - %.0f A/m²\n', J_min, J_max);
fprintf(fid, '  T_C:      %.1f - %.1f °C\n', T_min, T_max);
fprintf(fid, '  P:        %.2f - %.2f bar\n\n', P_min, P_max);

fprintf(fid, 'FIGURE 1 & 2: Pareto Front Objectives\n');
fprintf(fid, '  Energy:   %.4f - %.4f kWh/kg\n', Energy_min, Energy_max);
fprintf(fid, '  Cost:     %.4f - %.4f m²\n', Cost_min, Cost_max);
fprintf(fid, '  Stress:   %.4f - %.4f\n\n', Stress_min, Stress_max);

fprintf(fid, 'FIGURE 4: Correlation Matrices\n');
fprintf(fid, '  Energy-Cost:   %.3f\n', corr_objs(1,2));
fprintf(fid, '  Energy-Stress: %.3f\n', corr_objs(1,3));
fprintf(fid, '  Cost-Stress:   %.3f\n', corr_objs(2,3));
fprintf(fid, '  J-T:           %.3f\n\n', corr_decs(3,4));

fprintf(fid, 'TABLE 6: Representative Solutions\n');
fprintf(fid, '  Min Energy:  Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(minEnergyIdx,:));
fprintf(fid, '  Min Cost:    Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(minCostIdx,:));
fprintf(fid, '  Min Stress:  Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(minStressIdx,:));
fprintf(fid, '  Balanced:    Energy=%.2f, Cost=%.3f, Stress=%.2f\n', pareto_final(balancedIdx,:));

fclose(fid);

fprintf('\n✅ Values saved to: figure_values_for_captions.txt\n');
fprintf('%s\n', repmat('=', 1, 60));