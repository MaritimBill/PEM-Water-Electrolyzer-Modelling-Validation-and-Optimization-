%% ABLATION STUDY FOR FeaNSGA-II (CORRECTED)
clear; clc; close all;

%% ========== SETTINGS ==========
popSize = 100;
maxGen = 100;
maxFE = popSize * maxGen;
numRuns = 30;  % 30 runs for statistical significance
problemName = 'PEM_Problem_Scaled';

% Define variants with FLAGS
variants = {
    'All components',           true,  true,  true,  true,  true,  true;
    'No FeasibilityFirst',     false, true,  true,  true,  true,  true;
    'No EliteArchive',          true, false, true,  true,  true,  true;
    'No Immigrants',            true, true,  false, true,  true,  true;
    'No IntegerRepair',         true, true,  true,  false, true,  true;
    'Penalty Constraint',       true, true,  true,  true,  false, true;
    'Random Immigrants',        true, true,  'rand', true,  true,  true;
    'Aggressive Immigrants',    true, true,  0.20,  true,  true,  true;
    'Conservative Immigrants',  true, true,  0.02,  true,  true,  true;
    'Only FeasFirst',          false, false, false, true,  true,  false;
};

nVariants = size(variants, 1);
displayNames = variants(:,1);

% Add PlatEMO - UPDATE THIS PATH TO YOUR LOCATION
platemo_path = 'C:\Users\Admin\Desktop\PlatEMO-4.7';
if exist(platemo_path, 'dir')
    addpath(genpath(platemo_path));
else
    % Try current directory
    if exist('PlatEMO', 'dir')
        addpath(genpath('PlatEMO'));
    else
        error('PlatEMO not found! Please update platemo_path variable.');
    end
end

%% ========== STORAGE ==========
all_objs = {};
variant_objs = cell(nVariants, 1);
variant_times = cell(nVariants, 1);
variant_feas = cell(nVariants, 1);

fprintf('\n=== ABLATION STUDY: FeaNSGA-II ===\n');
fprintf('Variants: %d | Runs: %d | Generations: %d\n\n', nVariants, numRuns, maxGen);

%% ========== MAIN LOOP ==========
for run = 1:numRuns
    fprintf('Run %2d/%d: ', run, numRuns);
    rng(run*100, 'twister');
    
    for v = 1:nVariants
        use_feas_first = variants{v,2};
        use_archive = variants{v,3};
        immigrant_mode = variants{v,4};
        use_int_repair = variants{v,5};
        use_feas_selector = variants{v,6};
        use_archive_injection = variants{v,7};
        
        try
            Problem = feval(problemName);
            Problem.N = popSize;
            Problem.maxFE = maxFE;
            
            [objs, runtime, feas_count] = runFeaNSGA_II_with_flags(...
                Problem, popSize, maxGen, ...
                use_feas_first, use_archive, immigrant_mode, ...
                use_int_repair, use_feas_selector, use_archive_injection);
            
            variant_objs{v}{run} = objs;
            variant_times{v}(run) = runtime;
            variant_feas{v}(run) = feas_count / popSize * 100;
            
            if ~isempty(objs)
                all_objs{end+1} = objs;
                fprintf('.');
            else
                fprintf('x');
            end
            
        catch ME
            variant_objs{v}{run} = [];
            variant_times{v}(run) = NaN;
            variant_feas{v}(run) = 0;
            fprintf('E');
            if run == 1 && v == 1
                fprintf('\nError: %s\n', ME.message);
            end
        end
    end
    fprintf(' done\n');
end

%% ========== BUILD REFERENCE FRONT ==========
fprintf('\nBuilding reference front...\n');

all_objs_matrix = [];
for i = 1:length(all_objs)
    if ~isempty(all_objs{i})
        all_objs_matrix = [all_objs_matrix; all_objs{i}];
    end
end

if isempty(all_objs_matrix)
    error('No feasible solutions found!');
end

fprintf('Total solutions: %d\n', size(all_objs_matrix, 1));

ND = NDSort(all_objs_matrix, size(all_objs_matrix,1)) == 1;
refFront = all_objs_matrix(ND, :);
fprintf('Reference front size: %d points\n', size(refFront, 1));

minF = min(refFront, [], 1);
maxF = max(refFront, [], 1);
refFront_norm = (refFront - minF) ./ (maxF - minF + 1e-12);
refPoint = max(refFront_norm, [], 1) * 1.1;

%% ========== CALCULATE METRICS ==========
fprintf('Calculating metrics...\n');

hv_vals = zeros(nVariants, 1);
igd_vals = zeros(nVariants, 1);
time_vals = zeros(nVariants, 1);
feas_vals = zeros(nVariants, 1);
hv_std = zeros(nVariants, 1);
igd_std = zeros(nVariants, 1);
time_std = zeros(nVariants, 1);

% Store per-run HV for boxplots
hv_per_run = cell(nVariants, 1);

for v = 1:nVariants
    all_run_hv = [];
    all_run_igd = [];
    
    for run = 1:numRuns
        objs = variant_objs{v}{run};
        
        if isempty(objs) || size(objs, 1) < 2
            all_run_hv = [all_run_hv; 0];
            all_run_igd = [all_run_igd; 1e6];
            continue;
        end
        
        objs_n = (objs - minF) ./ (maxF - minF + 1e-12);
        
        % HV (Monte Carlo)
        Nsample = 5000;
        samples = rand(Nsample, 3) .* refPoint;
        dominated = false(Nsample, 1);
        for i = 1:size(objs_n, 1)
            dominated = dominated | all(samples >= objs_n(i, :), 2);
        end
        hv = sum(dominated) / Nsample * prod(refPoint);
        all_run_hv = [all_run_hv; hv];
        
        % IGD
        D = pdist2(refFront_norm, objs_n);
        igd = mean(min(D, [], 2));
        all_run_igd = [all_run_igd; igd];
    end
    
    hv_per_run{v} = all_run_hv;
    hv_vals(v) = mean(all_run_hv);
    hv_std(v) = std(all_run_hv);
    igd_vals(v) = mean(all_run_igd);
    igd_std(v) = std(all_run_igd);
    time_vals(v) = mean(variant_times{v});
    time_std(v) = std(variant_times{v});
    feas_vals(v) = mean(variant_feas{v});
end

%% ========== PRINT RESULTS TABLE ==========
fprintf('\n\n========================================\n');
fprintf('ABLATION STUDY RESULTS\n');
fprintf('========================================\n\n');

fprintf('%-35s %15s %15s %10s %10s\n', ...
    'Variant', 'HV (↑)', 'IGD (↓)', 'Time(s)', 'Feas%%');
fprintf('%s\n', repmat('-', 95, 1));

for v = 1:nVariants
    if hv_vals(v) == 0
        fprintf('%-35s %15s %15s %8.2f %9.1f%%\n', ...
            displayNames{v}, 'FAILED', 'FAILED', time_vals(v), feas_vals(v));
    else
        fprintf('%-35s %8.4e ±%-6.2e %8.4e ±%-6.2e %8.2f %9.1f%%\n', ...
            displayNames{v}, hv_vals(v), hv_std(v), igd_vals(v), igd_std(v), time_vals(v), feas_vals(v));
    end
end

%% ========== FIND BEST VARIANTS ==========
[best_hv, best_idx] = max(hv_vals);
fprintf('\n🏆 BEST BY HV: %s (HV = %.4e)\n', displayNames{best_idx}, best_hv);

[best_igd, best_idx_igd] = min(igd_vals);
fprintf('🏆 BEST BY IGD: %s (IGD = %.4e)\n', displayNames{best_idx_igd}, best_igd);

%% ========== FIGURE 1: HV BAR CHART WITH ERROR BARS ==========
fig1 = figure('Position', [50, 50, 1200, 600], 'Color', 'white');

% Sort by HV for better visualization
[~, sort_idx] = sort(hv_vals, 'descend');
hv_sorted = hv_vals(sort_idx);
hv_std_sorted = hv_std(sort_idx);
names_sorted = displayNames(sort_idx);

% Create bar chart
subplot(1,2,1);
bar_h = bar(hv_sorted);
bar_h.FaceColor = [0.2, 0.4, 0.7];
bar_h.EdgeColor = 'k';
bar_h.LineWidth = 1.2;
hold on;

% Add error bars
errorbar(1:nVariants, hv_sorted, hv_std_sorted, 'k.', 'LineWidth', 1.2, 'MarkerSize', 10);

% Add baseline reference line
baseline_idx = find(strcmp(names_sorted, 'All components'));
if ~isempty(baseline_idx)
    yline(hv_sorted(baseline_idx), 'r--', 'LineWidth', 1.5);
    text(nVariants-2, hv_sorted(baseline_idx)*1.01, 'Baseline', 'Color', 'r', 'FontSize', 10);
end

xlabel('Variant', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ablation Study: HV Comparison', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XTick', 1:nVariants, 'XTickLabel', names_sorted, 'XTickLabelRotation', 45);
grid on;

% Add value labels on bars
for i = 1:length(hv_sorted)
    if hv_sorted(i) > 0
        text(i, hv_sorted(i) + hv_std_sorted(i)*0.5, sprintf('%.4f', hv_sorted(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

%% ========== FIGURE 1B: IGD BAR CHART ==========
subplot(1,2,2);
[~, sort_idx_igd] = sort(igd_vals, 'ascend');
igd_sorted = igd_vals(sort_idx_igd);
igd_std_sorted = igd_std(sort_idx_igd);
names_sorted_igd = displayNames(sort_idx_igd);

bar_h2 = bar(igd_sorted);
bar_h2.FaceColor = [0.7, 0.3, 0.2];
bar_h2.EdgeColor = 'k';
bar_h2.LineWidth = 1.2;
hold on;

errorbar(1:nVariants, igd_sorted, igd_std_sorted, 'k.', 'LineWidth', 1.2, 'MarkerSize', 10);

xlabel('Variant', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Inverted Generational Distance (IGD)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ablation Study: IGD Comparison', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XTick', 1:nVariants, 'XTickLabel', names_sorted_igd, 'XTickLabelRotation', 45);
grid on;

sgtitle('FeaNSGA-II Ablation Study: Performance Metrics', 'FontSize', 16, 'FontWeight', 'bold');

%% ========== FIGURE 2: BOXPLOT OF HV DISTRIBUTIONS ==========
fig2 = figure('Position', [100, 100, 1000, 600], 'Color', 'white');

% Prepare data for boxplot (only variants with some HV values)
hv_data = [];
group_labels = {};
for v = 1:nVariants
    if max(hv_per_run{v}) > 0
        hv_data = [hv_data; hv_per_run{v}];
        group_labels = [group_labels; repmat(displayNames(v), length(hv_per_run{v}), 1)];
    end
end

% Create boxplot
boxplot(hv_data, group_labels, 'Orientation', 'horizontal', 'Colors', 'k', ...
    'Symbol', 'o', 'OutlierSize', 4);
xlabel('Hypervolume (HV)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Variant', 'FontSize', 12, 'FontWeight', 'bold');
title('HV Distribution Across 30 Runs', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Highlight baseline
ax = gca;
ax.FontSize = 9;

% Add vertical line for baseline median
baseline_median = median(hv_per_run{1});
xline(baseline_median, 'r--', 'LineWidth', 1.5);
text(baseline_median*1.01, nVariants-2, 'Baseline Median', 'Color', 'r', 'FontSize', 10);

%% ========== FIGURE 3: COMPONENT IMPACT HEATMAP ==========
fig3 = figure('Position', [150, 150, 800, 600], 'Color', 'white');

% Calculate impact matrix
components = {
    'Feasibility-first', 2;
    'Elite Archive', 3;
    'Periodic Immigrants', 4;
    'Integer Repair', 5;
    'Penalty Constraint', 6;
    'Random Immigrants', 7;
    'Aggressive (20%%)', 8;
    'Conservative (2%%)', 9;
};

nComp = length(components);
impact_matrix = zeros(nComp, 2);  % HV impact, IGD impact

baseline_hv = hv_vals(1);
baseline_igd = igd_vals(1);

for i = 1:nComp
    idx = components{i,2};
    if hv_vals(idx) > 0 && baseline_hv > 0
        impact_matrix(i,1) = (hv_vals(idx) - baseline_hv) / baseline_hv * 100;
        impact_matrix(i,2) = (igd_vals(idx) - baseline_igd) / baseline_igd * 100;
    end
end

% Create heatmap
imagesc(impact_matrix);
colormap(redbluecmap());
colorbar;
caxis([-10, 10]);

set(gca, 'XTick', 1:2, 'XTickLabel', {'HV Impact (%)', 'IGD Impact (%)'});
set(gca, 'YTick', 1:nComp, 'YTickLabel', components(:,1));
title('Component Impact Analysis', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Metric', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Component Modified/Removed', 'FontSize', 12, 'FontWeight', 'bold');

% Add value labels
for i = 1:nComp
    for j = 1:2
        if impact_matrix(i,j) ~= 0
            text(j, i, sprintf('%.1f%%', impact_matrix(i,j)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    end
end

%% ========== FIGURE 4: RUNTIME COMPARISON ==========
fig4 = figure('Position', [200, 200, 800, 500], 'Color', 'white');

[~, sort_idx_time] = sort(time_vals, 'ascend');
time_sorted = time_vals(sort_idx_time);
names_sorted_time = displayNames(sort_idx_time);

bar(time_sorted, 'FaceColor', [0.3, 0.6, 0.4], 'EdgeColor', 'k');
xlabel('Variant', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Runtime (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
title('Computational Time Comparison', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XTick', 1:nVariants, 'XTickLabel', names_sorted_time, 'XTickLabelRotation', 45);
grid on;

% Add value labels
for i = 1:length(time_sorted)
    if ~isnan(time_sorted(i))
        text(i, time_sorted(i) + 0.02, sprintf('%.2fs', time_sorted(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
end

%% ========== SAVE ALL FIGURES ==========
saveas(fig1, 'ablation_HV_IGD_bars.png');
saveas(fig2, 'ablation_HV_boxplot.png');
saveas(fig3, 'ablation_impact_heatmap.png');
saveas(fig4, 'ablation_runtime.png');

fprintf('\n✅ Figures saved:\n');
fprintf('   • ablation_HV_IGD_bars.png\n');
fprintf('   • ablation_HV_boxplot.png\n');
fprintf('   • ablation_impact_heatmap.png\n');
fprintf('   • ablation_runtime.png\n');

%% ========== SAVE RESULTS ==========
save('ablation_results.mat', 'hv_vals', 'igd_vals', 'time_vals', 'feas_vals', ...
     'hv_std', 'igd_std', 'time_std', 'displayNames', 'hv_per_run');
fprintf('✓ Results saved to ablation_results.mat\n');

%% ========== HELPER: RED-BLUE COLORMAP ==========
function cmap = redbluecmap()
    n = 64;
    top = [0.9, 0.2, 0.2];   % Red (negative impact)
    middle = [1, 1, 1];       % White (neutral)
    bottom = [0.2, 0.2, 0.9]; % Blue (positive impact)
    
    cmap = [linspace(top(1), middle(1), n/2)', ...
            linspace(top(2), middle(2), n/2)', ...
            linspace(top(3), middle(3), n/2)';
            linspace(middle(1), bottom(1), n/2)', ...
            linspace(middle(2), bottom(2), n/2)', ...
            linspace(middle(3), bottom(3), n/2)'];
end

%% ==================== ALGORITHM FUNCTIONS ====================

function [objs, runtime, feas_count] = runFeaNSGA_II_with_flags(...
    Problem, N, maxGen, use_feas_first, use_archive, immigrant_mode, ...
    use_int_repair, use_feas_selector, use_archive_injection)
    
    D = Problem.D;
    tic;
    
    Population = Problem.Initialization();
    Population = repairAndEvaluate(Problem, Population, use_int_repair);
    [Population, FrontNo, CrowdDis] = EnvironmentalSelection(Population, N);
    
    Archive = [];
    if use_archive
        Archive = getFeasibleNondominated(Population);
        if isempty(Archive), Archive = Population; end
    end
    
    gen = 0;
    while gen < maxGen
        CV = constraintViolation(Population);
        
        if use_feas_first
            MatingPool = zeros(1, N);
            for i = 1:N
                a = randi(N); b = randi(N);
                MatingPool(i) = tournamentWinner(a, b, FrontNo, CrowdDis, CV);
            end
        else
            MatingPool = TournamentSelection(2, N, FrontNo, CrowdDis);
        end
        Parents = Population(MatingPool);
        
        if use_archive && use_archive_injection && ~isempty(Archive)
            eliteNum = min(max(2, round(0.15*N)), numel(Archive));
            eliteIdx = randi(numel(Archive), 1, eliteNum);
            replaceIdx = randperm(N, eliteNum);
            Parents(replaceIdx) = Archive(eliteIdx);
        end
        
        Offspring = OperatorGA(Problem, Parents);
        Offspring = repairAndEvaluate(Problem, Offspring, use_int_repair);
        
        if ischar(immigrant_mode) && strcmp(immigrant_mode, 'rand')
            if rand < 0.3
                k = max(2, round(0.08*N));
                ImmDec = repmat(Problem.lower, k, 1) + rand(k, D) .* repmat(Problem.upper - Problem.lower, k, 1);
                if D >= 2 && use_int_repair, ImmDec(:,2) = round(ImmDec(:,2)); end
                ImmPop = Problem.Evaluation(ImmDec);
                Offspring = [Offspring, ImmPop];
            end
        elseif isnumeric(immigrant_mode) && immigrant_mode > 0
            if mod(gen, 10) == 0
                k = max(2, round(immigrant_mode * N));
                ImmDec = repmat(Problem.lower, k, 1) + rand(k, D) .* repmat(Problem.upper - Problem.lower, k, 1);
                if D >= 2 && use_int_repair, ImmDec(:,2) = round(ImmDec(:,2)); end
                ImmPop = Problem.Evaluation(ImmDec);
                Offspring = [Offspring, ImmPop];
            end
        elseif immigrant_mode == true
            if mod(gen, 10) == 0
                k = max(2, round(0.08*N));
                ImmDec = repmat(Problem.lower, k, 1) + rand(k, D) .* repmat(Problem.upper - Problem.lower, k, 1);
                if D >= 2 && use_int_repair, ImmDec(:,2) = round(ImmDec(:,2)); end
                ImmPop = Problem.Evaluation(ImmDec);
                Offspring = [Offspring, ImmPop];
            end
        end
        
        [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], N);
        
        if use_archive
            NewArchive = getFeasibleNondominated(Population);
            if ~isempty(NewArchive)
                Archive = mergeArchive(Archive, NewArchive, N);
            end
        end
        
        gen = gen + 1;
    end
    
    runtime = toc;
    objs = Population.objs;
    cons = Population.cons;
    if ~isempty(cons)
        feasible = all(cons <= 1e-6, 2);
        objs = objs(feasible, :);
        feas_count = sum(feasible);
    else
        feas_count = length(Population);
    end
end

function Pop = repairAndEvaluate(Problem, Pop, use_int_repair)
    if isempty(Pop), return; end
    Dec = Pop.decs;
    LB = Problem.lower;
    UB = Problem.upper;
    Dec = max(repmat(LB, size(Dec,1), 1), min(repmat(UB, size(Dec,1), 1), Dec));
    if use_int_repair && size(Dec,2) >= 2
        Dec(:,2) = round(Dec(:,2));
        Dec(:,2) = max(LB(2), min(UB(2), Dec(:,2)));
    end
    bad = any(~isfinite(Dec),2);
    if any(bad)
        m = sum(bad);
        Dec(bad,:) = repmat(LB, m, 1) + rand(m, Problem.D) .* repmat(UB-LB, m, 1);
        if use_int_repair && size(Dec,2) >= 2
            Dec(bad,2) = round(Dec(bad,2));
        end
    end
    Pop = Problem.Evaluation(Dec);
end

function CV = constraintViolation(Pop)
    if isempty(Pop), CV = []; return; end
    Con = Pop.cons;
    if isempty(Con), CV = zeros(numel(Pop),1);
    else, CV = sum(max(0,Con),2); end
end

function winner = tournamentWinner(a, b, FrontNo, CrowdDis, CV)
    feasA = CV(a) <= 1e-12;
    feasB = CV(b) <= 1e-12;
    if feasA && ~feasB, winner = a;
    elseif ~feasA && feasB, winner = b;
    elseif feasA && feasB
        if FrontNo(a) < FrontNo(b), winner = a;
        elseif FrontNo(b) < FrontNo(a), winner = b;
        else
            if CrowdDis(a) > CrowdDis(b), winner = a;
            else, winner = b; end
        end
    else
        if CV(a) < CV(b), winner = a;
        else, winner = b; end
    end
end

function Feas = isFeasible(Pop)
    if isempty(Pop), Feas = false(0,1); return; end
    Con = Pop.cons;
    if isempty(Con), Feas = true(numel(Pop),1);
    else, Feas = all(Con <= 1e-12, 2); end
end

function Archive = getFeasibleNondominated(Pop)
    Feas = isFeasible(Pop);
    Archive = Pop(Feas);
    if isempty(Archive), return; end
    [FrontNo,~] = NDSort(Archive.objs,1);
    Archive = Archive(FrontNo == 1);
end

function Archive = mergeArchive(A, B, N)
    Archive = [A, B];
    if isempty(Archive), return; end
    Feas = isFeasible(Archive);
    Archive = Archive(Feas);
    if isempty(Archive), return; end
    [FrontNo,~] = NDSort(Archive.objs,1);
    Archive = Archive(FrontNo == 1);
    cap = max(3, 5*N);
    if numel(Archive) > cap
        idx = randperm(numel(Archive), cap);
        Archive = Archive(idx);
    end
end