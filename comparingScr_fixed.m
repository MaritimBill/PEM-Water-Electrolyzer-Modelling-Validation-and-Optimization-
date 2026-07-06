%% FIXED COMPREHENSIVE 16-ALGORITHM COMPARISON (FAIR + CORRECT HV)
% With proper handling of infeasible solutions, feasibility tracking, and time
% Includes individual rankings for each metric
% ADDED: Mean ± std for all metrics
% ADDED: Wilcoxon signed-rank tests with Holm-Bonferroni correction

clear; clc; close all;

popSize = 100; maxGen = 100; maxFE = popSize * maxGen; numRuns = 30;  % Increased to 30 for statistical validity

algorithms = {'NSGAII','NSGAIII','ANSGAIII','GDE3','MOEAD','MOEADDE','MOEADSTM','ARMOEA','RVEA','RSEA','SRA','KnEA','BCEIBEA','SPEA2','IBEA','FeaNSGA_II'};
displayNames = {'NSGA-II','NSGA-III','ANSGA-III','GDE3','MOEA/D','MOEA/D-DE','MOEA/D-STM','AR-MOEA','RVEA','RSEA','SRA','KnEA','BCE-IBEA','SPEA2','IBEA','FeaNSGA_II ◊'};

nAlgos = length(algorithms);

% Store per-run metrics for statistical analysis
hv_per_run = cell(nAlgos, 1);
igd_per_run = cell(nAlgos, 1);
gd_per_run = cell(nAlgos, 1);
spread_per_run = cell(nAlgos, 1);
time_per_run = cell(nAlgos, 1);
feas_per_run = cell(nAlgos, 1);

results = struct;
all_feasible_objs = [];

% Track additional metrics
feasible_count = zeros(nAlgos, 1);      % Number of feasible solutions across runs
total_runs = zeros(nAlgos, 1);          % Number of successful runs
total_time = zeros(nAlgos, 1);          % Total computation time
feasibility_rate = zeros(nAlgos, 1);    % Percentage of feasible solutions

fprintf('\n=== RUNNING ===\n');

%% RUN
for run=1:numRuns
    fprintf('\n========== RUN %d / %d ==========\n', run, numRuns);
    for a=1:nAlgos
        algoName=algorithms{a};
        rng(run*100+a,'twister');

        Problem=PEM_Problem_Scaled();
        Problem.N=popSize; Problem.maxFE=maxFE;

        fprintf('%-20s', displayNames{a});

        try
            Algorithm=feval(algoName); 
            t_start = tic;
            Algorithm.Solve(Problem);
            runtime = toc(t_start);

            finalPop=[];
            if isprop(Algorithm,'result') && ~isempty(Algorithm.result)
                finalPop=Algorithm.result{end,2};
            end

            if isempty(finalPop)
                fprintf(' X (no population)\n');
                total_runs(a) = total_runs(a) + 1;
                continue;
            end

            objs = cat(1, finalPop.objs);
            cons = cat(1, finalPop.cons);
            feasible = all(cons<=0,2);
            n_feasible = sum(feasible);

            if n_feasible == 0
                fprintf(' no feasible\n');
                total_runs(a) = total_runs(a) + 1;
                continue;
            end

            objs = objs(feasible,:);
            all_feasible_objs = [all_feasible_objs; objs];

            % Store results
            if ~isfield(results, algoName)
                results.(algoName).objs = [];
                results.(algoName).time = [];
                results.(algoName).feasible_count = 0;
                results.(algoName).per_run_objs = {};
            end
            
            results.(algoName).objs = [results.(algoName).objs; objs];
            results.(algoName).time(end+1) = runtime;
            results.(algoName).feasible_count = results.(algoName).feasible_count + n_feasible;
            results.(algoName).per_run_objs{end+1} = objs;
            
            total_runs(a) = total_runs(a) + 1;
            feasible_count(a) = feasible_count(a) + n_feasible;
            total_time(a) = total_time(a) + runtime;

            fprintf(' ok (%d feasible)\n', n_feasible);

        catch err
            fprintf(' error: %s\n', err.message);
            total_runs(a) = total_runs(a) + 1;
        end
    end
end

%% Calculate feasibility rates
for a=1:nAlgos
    if total_runs(a) > 0
        avg_pop_per_run = 100;  % Expected population size
        expected_total = total_runs(a) * avg_pop_per_run;
        feasibility_rate(a) = feasible_count(a) / expected_total * 100;
    else
        feasibility_rate(a) = 0;
    end
end

%% =========================
%% REFERENCE FRONT (FIXED)
%% =========================
ND = NDSort(all_feasible_objs,size(all_feasible_objs,1))==1;
refFront = all_feasible_objs(ND,:);

% Normalize reference front
minF = min(refFront,[],1);
maxF = max(refFront,[],1);
refFront_norm = (refFront - minF) ./ (maxF - minF + 1e-12);

%% =========================
%% METRICS (with per-run storage for statistics)
%% =========================

hv_vals = inf(nAlgos,1);
igd_vals = inf(nAlgos,1);
gd_vals = inf(nAlgos,1);
spread_vals = inf(nAlgos,1);
avg_time = zeros(nAlgos,1);
std_time = zeros(nAlgos,1);

% Per-run metric storage
hv_per_run_vals = cell(nAlgos,1);
igd_per_run_vals = cell(nAlgos,1);
gd_per_run_vals = cell(nAlgos,1);
spread_per_run_vals = cell(nAlgos,1);
feas_per_run_vals = cell(nAlgos,1);

refPoint = ones(1,3)*1.1; % normalized reference

for a=1:nAlgos
    algoName=algorithms{a};

    if ~isfield(results,algoName) || isempty(results.(algoName).objs)
        hv_vals(a) = 0;
        igd_vals(a) = 1e6;  % Penalty: very large
        gd_vals(a) = 1e6;    % Penalty: very large
        spread_vals(a) = 1e6; % Penalty: very large
        avg_time(a) = total_time(a) / max(1, total_runs(a));
        std_time(a) = 0;
        hv_per_run_vals{a} = [];
        igd_per_run_vals{a} = [];
        gd_per_run_vals{a} = [];
        spread_per_run_vals{a} = [];
        feas_per_run_vals{a} = [];
        continue;
    end
    
    objs = results.(algoName).objs;
    times = results.(algoName).time;
    avg_time(a) = mean(times);
    std_time(a) = std(times);
    
    % Store feasibility per run
    feas_per_run = zeros(total_runs(a), 1);
    for r = 1:total_runs(a)
        run_objs = results.(algoName).per_run_objs{r};
        feas_per_run(r) = size(run_objs, 1);
    end
    feas_per_run_vals{a} = feas_per_run;

    % Per-run metrics
    hv_run = zeros(total_runs(a), 1);
    igd_run = zeros(total_runs(a), 1);
    gd_run = zeros(total_runs(a), 1);
    spread_run = zeros(total_runs(a), 1);
    
    for r = 1:total_runs(a)
        run_objs = results.(algoName).per_run_objs{r};
        
        % Normalize
        objs_n = (run_objs - minF) ./ (maxF - minF + 1e-12);
        
        %% HV (Monte Carlo approx)
        if size(objs_n,1) > 0
            Nsample = 10000;
            samples = rand(Nsample,3).*refPoint;
            dominated = false(Nsample,1);
            for i=1:size(objs_n,1)
                dominated = dominated | all(samples >= objs_n(i,:),2);
            end
            hv_run(r) = sum(dominated)/Nsample * prod(refPoint);
        else
            hv_run(r) = 0;
        end
        
        %% IGD
        D = pdist2(refFront_norm, objs_n);
        igd_run(r) = mean(min(D,[],2));
        
        %% GD
        D2 = pdist2(objs_n, refFront_norm);
        gd_run(r) = mean(min(D2,[],2));
        
        %% Spread
        if size(objs_n,1) > 1
            [~,idx]=sort(objs_n(:,1));
            objs_s=objs_n(idx,:);
            d = sqrt(sum(diff(objs_s).^2,2));
            spread_run(r) = std(d)/(mean(d)+1e-12);
        else
            spread_run(r) = 1e6;
        end
    end
    
    hv_per_run_vals{a} = hv_run;
    igd_per_run_vals{a} = igd_run;
    gd_per_run_vals{a} = gd_run;
    spread_per_run_vals{a} = spread_run;
    
    % Store means
    hv_vals(a) = mean(hv_run);
    igd_vals(a) = mean(igd_run);
    gd_vals(a) = mean(gd_run);
    spread_vals(a) = mean(spread_run);
end

%% =========================
%% STATISTICAL TESTS (Wilcoxon with Holm-Bonferroni)
%% =========================

fprintf('\n\n================ STATISTICAL ANALYSIS =================\n\n');

% Identify best algorithm by HV (primary metric)
[best_hv_val, best_idx] = max(hv_vals);
best_algo = algorithms{best_idx};
best_display = displayNames{best_idx};

fprintf('Best algorithm by HV: %s (%.4f)\n\n', best_display, best_hv_val);

% Perform pairwise Wilcoxon signed-rank tests for each metric
metrics_names = {'HV', 'IGD', 'GD', 'Spread', 'Time'};
metrics_vals = {hv_vals, igd_vals, gd_vals, spread_vals, avg_time};
metrics_per_run = {hv_per_run_vals, igd_per_run_vals, gd_per_run_vals, spread_per_run_vals, time_per_run};
better_direction = {'higher', 'lower', 'lower', 'lower', 'lower'};  % What indicates better performance

for m = 1:length(metrics_names)
    metric_name = metrics_names{m};
    metric_vals = metrics_vals{m};
    metric_per_run = metrics_per_run{m};
    better = better_direction{m};
    
    fprintf('\n=== %s Metric (Wilcoxon signed-rank test, Holm-Bonferroni corrected) ===\n', metric_name);
    fprintf('Best algorithm: ');
    
    % Find best algorithm for this metric
    if strcmp(better, 'higher')
        [best_val, best_metric_idx] = max(metric_vals);
    else
        [best_val, best_metric_idx] = min(metric_vals);
    end
    fprintf('%s (%.4f)\n\n', displayNames{best_metric_idx}, best_val);
    
    % Collect p-values for all pairwise comparisons with best
    p_values = [];
    algo_indices = [];
    
    for a = 1:nAlgos
        if a == best_metric_idx || isempty(metric_per_run{a}) || isempty(metric_per_run{best_metric_idx})
            continue;
        end
        
        % Ensure same number of runs for comparison
        n_runs = min(length(metric_per_run{a}), length(metric_per_run{best_metric_idx}));
        if n_runs < 3
            continue;
        end
        
        % Get paired data
        data_best = metric_per_run{best_metric_idx}(1:n_runs);
        data_a = metric_per_run{a}(1:n_runs);
        
        % Wilcoxon signed-rank test
        try
            if strcmp(better, 'higher')
                % Best should have higher values
                p = signrank(data_a, data_best);  % Tests if data_a < data_best
            else
                % Best should have lower values
                p = signrank(data_best, data_a);  % Tests if data_best < data_a
            end
            p_values = [p_values; p];
            algo_indices = [algo_indices; a];
        catch
            % Skip if test fails
            continue;
        end
    end
    
    % Apply Holm-Bonferroni correction
    if ~isempty(p_values)
        [sorted_p, sort_idx] = sort(p_values);
        n_tests = length(sorted_p);
        alpha = 0.05;
        
        fprintf('%-20s %12s %12s %s\n', 'Algorithm', 'p-value', 'Adj. Threshold', 'Significant');
        fprintf('%s\n', repmat('-', 60, 1));
        
        for k = 1:n_tests
            adj_threshold = alpha / (n_tests - k + 1);
            is_significant = sorted_p(k) < adj_threshold;
            
            a_idx = algo_indices(sort_idx(k));
if is_significant; sig_str = '✓ YES'; else; sig_str = '✗ NO'; end            
            fprintf('%-20s %12.4e %12.4f %s\n', ...
                displayNames{a_idx}, sorted_p(k), adj_threshold, sig_str);
        end
    else
        fprintf('No valid comparisons available.\n');
    end
end

%% =========================
%% PRINT MAIN RESULTS (with Mean ± Std)
%% =========================

fprintf('\n\n================ FINAL =================\n');
fprintf('%-15s %12s %12s %12s %12s %10s %10s %s\n',...
    'Algorithm', 'HV', 'IGD', 'GD', 'Spr', 'Time(s)', 'Feas%', '');

for i=1:nAlgos
    a = i;  % Print in original order
    
    % Format feasibility display
    if feasibility_rate(a) == 0
        feas_str = '0.0% ❌';
    elseif feasibility_rate(a) < 50
        feas_str = sprintf('%.1f%% ⚠️', feasibility_rate(a));
    else
        feas_str = sprintf('%.1f%% ✓', feasibility_rate(a));
    end
    
    % Format with mean ± std
    if ~isempty(hv_per_run_vals{a}) && length(hv_per_run_vals{a}) > 0
        hv_std = std(hv_per_run_vals{a});
        igd_std = std(igd_per_run_vals{a});
        gd_std = std(gd_per_run_vals{a});
        spread_std = std(spread_per_run_vals{a});
        
        hv_str = sprintf('%.3f±%.3f', hv_vals(a), hv_std);
        igd_str = sprintf('%.3f±%.3f', igd_vals(a), igd_std);
        gd_str = sprintf('%.3f±%.3f', gd_vals(a), gd_std);
        spread_str = sprintf('%.3f±%.3f', spread_vals(a), spread_std);
        time_str = sprintf('%.2f±%.2f', avg_time(a), std_time(a));
    else
        hv_str = sprintf('%.3f', hv_vals(a));
        igd_str = sprintf('%.3f', igd_vals(a));
        gd_str = sprintf('%.3f', gd_vals(a));
        spread_str = sprintf('%.3f', spread_vals(a));
        time_str = sprintf('%.2f', avg_time(a));
    end
    
    % Mark if algorithm failed
    if hv_vals(a) == 0 && feasibility_rate(a) < 10
        fprintf('%-15s %12s %12s %12s %12s %10s %10s ❌ FAILED\n',...
            displayNames{a}, hv_str, igd_str, gd_str, spread_str, time_str, feas_str);
    else
        fprintf('%-15s %12s %12s %12s %12s %10s %10s\n',...
            displayNames{a}, hv_str, igd_str, gd_str, spread_str, time_str, feas_str);
    end
end

%% =========================
%% RANKING (WITH FEASIBILITY PENALTY)
%% =========================

% Rank each metric (lower rank = better)
[~, hv_rank] = sort(hv_vals, 'descend');
[~, igd_rank] = sort(igd_vals, 'ascend');
[~, gd_rank] = sort(gd_vals, 'ascend');
[~, spread_rank] = sort(spread_vals, 'ascend');
[~, time_rank] = sort(avg_time, 'ascend');
[~, feas_rank] = sort(feasibility_rate, 'descend');  % Higher feasibility = better rank

% Convert to rank numbers
rankHV = zeros(nAlgos,1);
rankIGD = zeros(nAlgos,1);
rankGD = zeros(nAlgos,1);
rankSpread = zeros(nAlgos,1);
rankTime = zeros(nAlgos,1);
rankFeas = zeros(nAlgos,1);

for i=1:nAlgos
    rankHV(hv_rank(i)) = i;
    rankIGD(igd_rank(i)) = i;
    rankGD(gd_rank(i)) = i;
    rankSpread(spread_rank(i)) = i;
    rankTime(time_rank(i)) = i;
    rankFeas(feas_rank(i)) = i;
end

% Weighted overall rank (with feasibility penalty for failed algorithms)
% Heavily penalize algorithms with low feasibility
feas_penalty = 1 ./ (feasibility_rate/100 + 0.01);  % 1 for 100% feasible, 100 for 0%

overall_rank = 0.35*rankHV + 0.20*rankIGD + 0.10*rankGD + ...
               0.10*rankSpread + 0.10*rankTime + 0.15*rankFeas;
overall_rank = overall_rank .* feas_penalty';  % Apply penalty

[~,order] = sort(overall_rank);

%% =========================
%% PRINT RANKINGS (with stats)
%% =========================

fprintf('\n\n================ INDIVIDUAL RANKINGS =================\n\n');

% HV Ranking
fprintf('=== HV RANKING (Higher is better) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'HV (Mean±Std)', 'Rank');
for i = 1:nAlgos
    a = hv_rank(i);
    if ~isempty(hv_per_run_vals{a}) && length(hv_per_run_vals{a}) > 0
        hv_str = sprintf('%.3f±%.3f', hv_vals(a), std(hv_per_run_vals{a}));
    else
        hv_str = sprintf('%.3f', hv_vals(a));
    end
    fprintf('%-20s %12s %6d\n', displayNames{a}, hv_str, i);
end

fprintf('\n=== IGD RANKING (Lower is better) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'IGD (Mean±Std)', 'Rank');
for i = 1:nAlgos
    a = igd_rank(i);
    if igd_vals(a) >= 1e5
        fprintf('%-20s %12s %6d ❌\n', displayNames{a}, '1e6', i);
    elseif ~isempty(igd_per_run_vals{a}) && length(igd_per_run_vals{a}) > 0
        igd_str = sprintf('%.3f±%.3f', igd_vals(a), std(igd_per_run_vals{a}));
        fprintf('%-20s %12s %6d\n', displayNames{a}, igd_str, i);
    else
        fprintf('%-20s %12.3f %6d\n', displayNames{a}, igd_vals(a), i);
    end
end

fprintf('\n=== GD RANKING (Lower is better) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'GD (Mean±Std)', 'Rank');
for i = 1:nAlgos
    a = gd_rank(i);
    if gd_vals(a) >= 1e5
        fprintf('%-20s %12s %6d ❌\n', displayNames{a}, '1e6', i);
    elseif ~isempty(gd_per_run_vals{a}) && length(gd_per_run_vals{a}) > 0
        gd_str = sprintf('%.3f±%.3f', gd_vals(a), std(gd_per_run_vals{a}));
        fprintf('%-20s %12s %6d\n', displayNames{a}, gd_str, i);
    else
        fprintf('%-20s %12.3f %6d\n', displayNames{a}, gd_vals(a), i);
    end
end

fprintf('\n=== SPREAD RANKING (Lower is better) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'Spread (Mean±Std)', 'Rank');
for i = 1:nAlgos
    a = spread_rank(i);
    if spread_vals(a) >= 1e5
        fprintf('%-20s %12s %6d ❌\n', displayNames{a}, '1e6', i);
    elseif ~isempty(spread_per_run_vals{a}) && length(spread_per_run_vals{a}) > 0
        spread_str = sprintf('%.3f±%.3f', spread_vals(a), std(spread_per_run_vals{a}));
        fprintf('%-20s %12s %6d\n', displayNames{a}, spread_str, i);
    else
        fprintf('%-20s %12.3f %6d\n', displayNames{a}, spread_vals(a), i);
    end
end

fprintf('\n=== TIME RANKING (Lower is better) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'Time (s) (Mean±Std)', 'Rank');
for i = 1:nAlgos
    a = time_rank(i);
    time_str = sprintf('%.2f±%.2f', avg_time(a), std_time(a));
    fprintf('%-20s %12s %6d\n', displayNames{a}, time_str, i);
end

fprintf('\n=== FEASIBILITY RANKING (Higher is better) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'Feasibility %', 'Rank');
for i = 1:nAlgos
    a = feas_rank(i);
    if feasibility_rate(a) == 0
        fprintf('%-20s %12.1f%% %6d ❌\n', displayNames{a}, feasibility_rate(a), i);
    elseif feasibility_rate(a) < 50
        fprintf('%-20s %12.1f%% %6d ⚠️\n', displayNames{a}, feasibility_rate(a), i);
    else
        fprintf('%-20s %12.1f%% %6d ✓\n', displayNames{a}, feasibility_rate(a), i);
    end
end

%% =========================
%% SUMMARY STATISTICS
%% =========================

fprintf('\n\n================ SUMMARY STATISTICS =================\n');
fprintf('Total runs per algorithm: %d\n', numRuns);
fprintf('Reference front size: %d points\n', size(refFront,1));

% Calculate and display top 3 algorithms by overall rank
fprintf('\n=== TOP 3 ALGORITHMS ===\n');
for i = 1:min(3, nAlgos)
    a = order(i);
    fprintf('%d. %s (Overall Score: %.2f)\n', i, displayNames{a}, overall_rank(a));
    
    % Format with std
    if ~isempty(hv_per_run_vals{a}) && length(hv_per_run_vals{a}) > 0
        fprintf('   HV: %.3f±%.3f | IGD: %.3f±%.3f | GD: %.3f±%.3f | Spr: %.3f±%.3f | Time: %.2f±%.2fs | Feas: %.1f%%\n', ...
            hv_vals(a), std(hv_per_run_vals{a}), ...
            igd_vals(a), std(igd_per_run_vals{a}), ...
            gd_vals(a), std(gd_per_run_vals{a}), ...
            spread_vals(a), std(spread_per_run_vals{a}), ...
            avg_time(a), std_time(a), feasibility_rate(a));
    else
        fprintf('   HV: %.3f | IGD: %.3f | GD: %.3f | Spr: %.3f | Time: %.2fs | Feas: %.1f%%\n', ...
            hv_vals(a), igd_vals(a), gd_vals(a), spread_vals(a), avg_time(a), feasibility_rate(a));
    end
end

% Identify failed algorithms
failed_algos = find(feasibility_rate == 0);
if ~isempty(failed_algos)
    fprintf('\n=== FAILED ALGORITHMS (0%% Feasibility) ===\n');
    for i = 1:length(failed_algos)
        a = failed_algos(i);
        fprintf('• %s\n', displayNames{a});
    end
end

fprintf('\n✅ Analysis complete!\n');

% Display statistical test summary
fprintf('\n=== STATISTICAL TEST SUMMARY ===\n');
fprintf('• Wilcoxon signed-rank tests performed for each metric\n');
fprintf('• Holm-Bonferroni correction applied to control family-wise error rate\n');
fprintf('• Significance level: α = 0.05\n');
fprintf('• Best algorithm by HV: %s (%.4f)\n', best_display, best_hv_val);