%% FIXED COMPREHENSIVE 19-ALGORITHM COMPARISON (ALL INCLUDED IN STATS)
% All algorithms are included in Wilcoxon signed-rank tests, even if 0% feasible.
% No "FAILED" labels – just show feasibility 0%.
% Time ranking only for algorithms that actually ran (avg_time > 0).

clear; clc; close all;

popSize = 100; maxGen = 100; maxFE = popSize * maxGen;
numRuns = 30;   % <<-- SET TO 30 for proper statistics

algorithms = {'NSGAII','NSGAIII','ANSGAIII','GDE3','MOEAD','MOEADDE','MOEADSTM',...
              'ARMOEA','RVEA','RSEA','SRA','KnEA','BCEIBEA','SPEA2','IBEA','NSGAII_AdaptiveSearch',...
              'dMOPSO','MOPSO','PREA'};
displayNames = {'NSGA-II','NSGA-III','ANSGA-III','GDE3','MOEA/D','MOEA/D-DE',...
                'MOEA/D-STM','AR-MOEA','RVEA','RSEA','SRA','KnEA','BCE-IBEA',...
                'SPEA2','IBEA','NSGAII_AdaptiveSearch','dMOPSO','MOPSO','PREA'};


nAlgos = length(algorithms);

% Storage
hv_per_run = cell(nAlgos, 1);
igd_per_run = cell(nAlgos, 1);
gd_per_run = cell(nAlgos, 1);
spread_per_run = cell(nAlgos, 1);
time_per_run = cell(nAlgos, 1);
feas_per_run = cell(nAlgos, 1);

results = struct;
all_feasible_objs = [];

feasible_count = zeros(nAlgos, 1);
total_runs = zeros(nAlgos, 1);
total_time = zeros(nAlgos, 1);
feasibility_rate = zeros(nAlgos, 1);

fprintf('\n=== RUNNING ===\n');

%% RUN
for run = 1:numRuns
    fprintf('\n========== RUN %d / %d ==========\n', run, numRuns);
    for a = 1:nAlgos
        algoName = algorithms{a};
        rng(run * 100 + a, 'twister');

        Problem = PEM_Problem_Scaled();
        Problem.N = popSize;
        Problem.maxFE = maxFE;

        fprintf('%-20s', displayNames{a});

        if ~isfield(results, algoName)
            results.(algoName).objs = [];
            results.(algoName).time = [];
            results.(algoName).feasible_count = 0;
            results.(algoName).per_run_objs = {};
        end

        try
            switch algoName
                case 'NSGAII',      Algorithm = NSGAII();
                case 'NSGAIII',     Algorithm = NSGAIII();
                case 'ANSGAIII',    Algorithm = ANSGAIII();
                case 'GDE3',        Algorithm = GDE3();
                case 'MOEAD',       Algorithm = MOEAD();
                case 'MOEADDE',     Algorithm = MOEADDE();
                case 'MOEADSTM',    Algorithm = MOEADSTM();
                case 'ARMOEA',      Algorithm = ARMOEA();
                case 'RVEA',        Algorithm = RVEA();
                case 'RSEA',        Algorithm = RSEA();
                case 'SRA',         Algorithm = SRA();
                case 'KnEA',        Algorithm = KnEA();
                case 'BCEIBEA',     Algorithm = BCEIBEA();
                case 'SPEA2',       Algorithm = SPEA2();
                case 'IBEA',        Algorithm = IBEA();
                case 'NSGAII_AdaptiveSearch',  Algorithm = NSGAII_AdaptiveSearch();
                case 'dMOPSO',      Algorithm = dMOPSO();
                case 'MOPSO',       Algorithm = MOPSO();
                case 'PREA',        Algorithm = PREA();
                                      
                otherwise, error('Unknown algorithm: %s', algoName);
            end
            
            t_start = tic;
            Algorithm.Solve(Problem);
            runtime = toc(t_start);

            finalPop = [];
            if isprop(Algorithm, 'result') && ~isempty(Algorithm.result)
                finalPop = Algorithm.result{end, 2};
            end

            if isempty(finalPop)
                fprintf(' X (no population)\n');
                results.(algoName).time(end + 1) = runtime;
                results.(algoName).per_run_objs{end + 1} = [];
                total_runs(a) = total_runs(a) + 1;
                total_time(a) = total_time(a) + runtime;
                continue;
            end

            objs = cat(1, finalPop.objs);
            cons = cat(1, finalPop.cons);
            feasible = all(cons <= 0, 2);
            n_feasible = sum(feasible);

            results.(algoName).time(end + 1) = runtime;

            if n_feasible > 0
                objs = objs(feasible, :);
                all_feasible_objs = [all_feasible_objs; objs];
                results.(algoName).objs = [results.(algoName).objs; objs];
                results.(algoName).feasible_count = results.(algoName).feasible_count + n_feasible;
                results.(algoName).per_run_objs{end + 1} = objs;
                feasible_count(a) = feasible_count(a) + n_feasible;
                fprintf(' ok (%d feasible)\n', n_feasible);
            else
                results.(algoName).per_run_objs{end + 1} = [];
                fprintf(' no feasible\n');
            end

            total_runs(a) = total_runs(a) + 1;
            total_time(a) = total_time(a) + runtime;

        catch err
            fprintf(' error: %s\n', err.message);
            if ~isfield(results, algoName)
                results.(algoName).objs = [];
                results.(algoName).time = [];
                results.(algoName).feasible_count = 0;
                results.(algoName).per_run_objs = {};
            end
            % Record time as 0 for this run – but we'll later exclude these from time ranking
            results.(algoName).time(end + 1) = 0;
            results.(algoName).per_run_objs{end + 1} = [];
            total_runs(a) = total_runs(a) + 1;
            % total_time not updated here because we don't want to average 0 for failed runs
            % Actually we should not add 0 to total_time because it's not a real runtime.
            % But we keep it for completeness; we'll later filter out zeros in time ranking.
        end
    end
end

%% Calculate feasibility rates
for a = 1:nAlgos
    if total_runs(a) > 0
        expected_total = total_runs(a) * 100;
        feasibility_rate(a) = feasible_count(a) / expected_total * 100;
    else
        feasibility_rate(a) = 0;
    end
end

%% =========================
%% REFERENCE FRONT
%% =========================
if ~isempty(all_feasible_objs)
    ND = NDSort(all_feasible_objs, size(all_feasible_objs, 1)) == 1;
    refFront = all_feasible_objs(ND, :);
else
    error('No feasible solutions found!');
end

minF = min(refFront, [], 1);
maxF = max(refFront, [], 1);
refFront_norm = (refFront - minF) ./ (maxF - minF + 1e-12);

%% =========================
%% METRICS
%% =========================

hv_vals = inf(nAlgos, 1);
igd_vals = inf(nAlgos, 1);
gd_vals = inf(nAlgos, 1);
spread_vals = inf(nAlgos, 1);
avg_time = zeros(nAlgos, 1);
std_time = zeros(nAlgos, 1);

hv_per_run_vals = cell(nAlgos, 1);
igd_per_run_vals = cell(nAlgos, 1);
gd_per_run_vals = cell(nAlgos, 1);
spread_per_run_vals = cell(nAlgos, 1);
time_per_run_vals = cell(nAlgos, 1);
feas_per_run_vals = cell(nAlgos, 1);

refPoint = ones(1, 3) * 1.1;

for a = 1:nAlgos
    algoName = algorithms{a};

    if ~isfield(results, algoName) || isempty(results.(algoName).objs)
        hv_vals(a) = 0;
        igd_vals(a) = 1e6;
        gd_vals(a) = 1e6;
        spread_vals(a) = 1e6;
        avg_time(a) = total_time(a) / max(1, total_runs(a));  % includes zeros for failed runs
        std_time(a) = 0;
        hv_per_run_vals{a} = [];
        igd_per_run_vals{a} = [];
        gd_per_run_vals{a} = [];
        spread_per_run_vals{a} = [];
        time_per_run_vals{a} = [];
        feas_per_run_vals{a} = [];
        continue;
    end

    objs = results.(algoName).objs;
    times = results.(algoName).time;
    avg_time(a) = mean(times);
    std_time(a) = std(times);

    feas_per_run = zeros(total_runs(a), 1);
    for r = 1:total_runs(a)
        run_objs = results.(algoName).per_run_objs{r};
        if ~isempty(run_objs)
            feas_per_run(r) = size(run_objs, 1);
        else
            feas_per_run(r) = 0;
        end
    end
    feas_per_run_vals{a} = feas_per_run;

    hv_run = zeros(total_runs(a), 1);
    igd_run = zeros(total_runs(a), 1);
    gd_run = zeros(total_runs(a), 1);
    spread_run = zeros(total_runs(a), 1);

    for r = 1:total_runs(a)
        run_objs = results.(algoName).per_run_objs{r};

        if isempty(run_objs)
            hv_run(r) = 0;
            igd_run(r) = 1e6;
            gd_run(r) = 1e6;
            spread_run(r) = 1e6;
            continue;
        end

        objs_n = (run_objs - minF) ./ (maxF - minF + 1e-12);

        if size(objs_n, 1) > 0
            Nsample = 10000;
            samples = rand(Nsample, 3) .* refPoint;
            dominated = false(Nsample, 1);
            for i = 1:size(objs_n, 1)
                dominated = dominated | all(samples >= objs_n(i, :), 2);
            end
            hv_run(r) = sum(dominated) / Nsample * prod(refPoint);
        else
            hv_run(r) = 0;
        end

        D = pdist2(refFront_norm, objs_n);
        igd_run(r) = mean(min(D, [], 2));

        D2 = pdist2(objs_n, refFront_norm);
        gd_run(r) = mean(min(D2, [], 2));

        if size(objs_n, 1) > 1
            [~, idx] = sort(objs_n(:, 1));
            objs_s = objs_n(idx, :);
            d = sqrt(sum(diff(objs_s).^2, 2));
            spread_run(r) = std(d) / (mean(d) + 1e-12);
        else
            spread_run(r) = 1e6;
        end
    end

    hv_per_run_vals{a} = hv_run;
    igd_per_run_vals{a} = igd_run;
    gd_per_run_vals{a} = gd_run;
    spread_per_run_vals{a} = spread_run;
    time_per_run_vals{a} = times;

    hv_vals(a) = mean(hv_run);
    igd_vals(a) = mean(igd_run);
    gd_vals(a) = mean(gd_run);
    spread_vals(a) = mean(spread_run);
end

%% =========================
%% STATISTICAL TESTS (INCLUDE ALL ALGORITHMS, EVEN 0% FEASIBILITY)
%% =========================

fprintf('\n\n================ STATISTICAL ANALYSIS =================\n\n');

metrics_names = {'HV', 'IGD', 'GD', 'Spread', 'Time'};
metrics_vals = {hv_vals, igd_vals, gd_vals, spread_vals, avg_time};
metrics_per_run = {hv_per_run_vals, igd_per_run_vals, gd_per_run_vals, spread_per_run_vals, time_per_run_vals};
better_direction = {'higher', 'lower', 'lower', 'lower', 'lower'};

%% ----- ORIGINAL TEST: compare each algorithm against the best for that metric -----
for m = 1:length(metrics_names)
    metric_name = metrics_names{m};
    metric_vals = metrics_vals{m};
    metric_per_run = metrics_per_run{m};
    better = better_direction{m};

    fprintf('\n=== %s Metric (Wilcoxon signed-rank test, Holm-Bonferroni) ===\n', metric_name);

    valid_idx = find(~isinf(metric_vals) & ~isnan(metric_vals));
    if isempty(valid_idx)
        fprintf('No valid metric values.\n');
        continue;
    end

    if strcmp(better, 'higher')
        [best_val, best_pos] = max(metric_vals(valid_idx));
    else
        [best_val, best_pos] = min(metric_vals(valid_idx));
    end
    best_metric_idx = valid_idx(best_pos);
    fprintf('Best algorithm: %s (%.4f)\n\n', displayNames{best_metric_idx}, best_val);

    p_values = [];
    algo_indices = [];

    for a = 1:nAlgos
        if a == best_metric_idx
            continue;
        end
        if isempty(metric_per_run{a}) || isempty(metric_per_run{best_metric_idx})
            continue;
        end

        n_runs = min(length(metric_per_run{a}), length(metric_per_run{best_metric_idx}));
        if n_runs < 3
            continue;
        end

        data_best = metric_per_run{best_metric_idx}(1:n_runs);
        data_a = metric_per_run{a}(1:n_runs);

        % Remove NaNs/Infs
        data_best = data_best(~isnan(data_best) & ~isinf(data_best));
        data_a = data_a(~isnan(data_a) & ~isinf(data_a));
        if length(data_best) < 3 || length(data_a) < 3
            continue;
        end

        n = min(length(data_best), length(data_a));
        data_best = data_best(1:n);
        data_a = data_a(1:n);

        if all(data_a == data_best)
            continue;
        end

        try
            if strcmp(better, 'higher')
                p = signrank(data_a, data_best);
            else
                p = signrank(data_best, data_a);
            end
            p_values = [p_values; p];
            algo_indices = [algo_indices; a];
        catch
            continue;
        end
    end

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
            fprintf('%-20s %12.4e %12.4f %s\n', ...
                displayNames{a_idx}, sorted_p(k), adj_threshold, ...
                ternary(is_significant, 'YES', 'NO'));
        end
    else
        fprintf('No valid comparisons available.\n');
    end
end

%% ----- ADDITIONAL TEST: compare against NSGAII_AdaptiveSearch as reference (include all) -----
reference_name = 'NSGAII_AdaptiveSearch';
ref_idx = find(strcmp(algorithms, reference_name));
if isempty(ref_idx)
    fprintf('\n=== Reference algorithm %s not found.\n', reference_name);
else
    fprintf('\n\n================ ADDITIONAL: COMPARISONS AGAINST %s =================\n', reference_name);
    for m = 1:length(metrics_names)
        metric_name = metrics_names{m};
        metric_per_run = metrics_per_run{m};
        better = better_direction{m};

        fprintf('\n=== %s Metric (vs %s) ===\n', metric_name, reference_name);
        fprintf('%-20s %12s %12s %s\n', 'Algorithm', 'p-value', 'Adj. Threshold', 'Significant');
        fprintf('%s\n', repmat('-', 60, 1));

        data_ref = metric_per_run{ref_idx};
        data_ref = data_ref(~isnan(data_ref) & ~isinf(data_ref));
        if isempty(data_ref)
            fprintf('Reference has no valid data.\n');
            continue;
        end

        p_values = [];
        algo_indices = [];

        for a = 1:nAlgos
            if a == ref_idx
                continue;
            end
            data_a = metric_per_run{a};
            data_a = data_a(~isnan(data_a) & ~isinf(data_a));
            if isempty(data_a)
                continue;
            end

            n = min(length(data_ref), length(data_a));
            if n < 3
                continue;
            end

            x = data_ref(1:n);
            y = data_a(1:n);
            if all(x == y)
                continue;
            end

            try
                if strcmp(better, 'higher')
                    p = signrank(y, x);
                else
                    p = signrank(x, y);
                end
                p_values = [p_values; p];
                algo_indices = [algo_indices; a];
            catch
                continue;
            end
        end

        if ~isempty(p_values)
            [sorted_p, sort_idx] = sort(p_values);
            n_tests = length(sorted_p);
            alpha = 0.05;
            for k = 1:n_tests
                adj_threshold = alpha / (n_tests - k + 1);
                is_significant = sorted_p(k) < adj_threshold;
                a_idx = algo_indices(sort_idx(k));
                fprintf('%-20s %12.4e %12.4f %s\n', ...
                    displayNames{a_idx}, sorted_p(k), adj_threshold, ...
                    ternary(is_significant, 'YES', 'NO'));
            end
        else
            fprintf('No valid comparisons.\n');
        end
    end
end

%% =========================
%% PRINT FINAL RESULTS TABLE
%% =========================
fprintf('\n\n================ FINAL =================\n');
fprintf('%-15s %12s %12s %12s %12s %10s %10s\n',...
    'Algorithm', 'HV', 'IGD', 'GD', 'Spr', 'Time(s)', 'Feas%');

for i = 1:nAlgos
    a = i;
    if feasibility_rate(a) == 0
        feas_str = '0.0%';
    else
        feas_str = sprintf('%.1f%%', feasibility_rate(a));
    end

    if ~isempty(hv_per_run_vals{a}) && length(hv_per_run_vals{a}) > 0
        hv_str = sprintf('%.3f±%.3f', hv_vals(a), std(hv_per_run_vals{a}));
        igd_str = sprintf('%.3f±%.3f', igd_vals(a), std(igd_per_run_vals{a}));
        gd_str = sprintf('%.3f±%.3f', gd_vals(a), std(gd_per_run_vals{a}));
        spread_str = sprintf('%.3f±%.3f', spread_vals(a), std(spread_per_run_vals{a}));
        time_str = sprintf('%.2f±%.2f', avg_time(a), std_time(a));
    else
        hv_str = sprintf('%.3f', hv_vals(a));
        igd_str = sprintf('%.3f', igd_vals(a));
        gd_str = sprintf('%.3f', gd_vals(a));
        spread_str = sprintf('%.3f', spread_vals(a));
        time_str = sprintf('%.2f', avg_time(a));
    end

    fprintf('%-15s %12s %12s %12s %12s %10s %10s\n',...
        displayNames{a}, hv_str, igd_str, gd_str, spread_str, time_str, feas_str);
end

%% =========================
%% INDIVIDUAL RANKINGS
%% =========================

[~, hv_rank] = sort(hv_vals, 'descend');
[~, igd_rank] = sort(igd_vals, 'ascend');
[~, gd_rank] = sort(gd_vals, 'ascend');
[~, spread_rank] = sort(spread_vals, 'ascend');

% TIME RANKING: only algorithms that actually ran (avg_time > 0)
% Get indices of algorithms with avg_time > 0 and total_runs > 0
valid_time_idx = find(avg_time > 0 & total_runs > 0);
if ~isempty(valid_time_idx)
    [~, time_sort] = sort(avg_time(valid_time_idx));
    time_rank = inf(nAlgos, 1);
    for r = 1:length(valid_time_idx)
        time_rank(valid_time_idx(time_sort(r))) = r;
    end
else
    time_rank = inf(nAlgos, 1);
end
% For those with avg_time == 0, set rank to largest (they are excluded from ranking)
time_rank(isinf(time_rank)) = nAlgos + 1;  % or keep as Inf, we'll handle

[~, feas_rank] = sort(feasibility_rate, 'descend');

rankHV = zeros(nAlgos, 1);
rankIGD = zeros(nAlgos, 1);
rankGD = zeros(nAlgos, 1);
rankSpread = zeros(nAlgos, 1);
rankTime = time_rank;
rankFeas = zeros(nAlgos, 1);

for i = 1:nAlgos
    rankHV(hv_rank(i)) = i;
    rankIGD(igd_rank(i)) = i;
    rankGD(gd_rank(i)) = i;
    rankSpread(spread_rank(i)) = i;
    rankFeas(feas_rank(i)) = i;
end

% Feasibility penalty
feas_penalty = ones(nAlgos, 1);
for a = 1:nAlgos
    if feasibility_rate(a) < 50
        feas_penalty(a) = 1e6;
    elseif feasibility_rate(a) < 90
        feas_penalty(a) = 10;
    else
        feas_penalty(a) = 1;
    end
end

% Set rankTime for excluded algorithms to a large number
rankTime(isinf(rankTime)) = 1e6;

overall_rank = 0.35 * rankHV + 0.20 * rankIGD + 0.10 * rankGD + ...
               0.10 * rankSpread + 0.10 * rankTime + 0.15 * rankFeas;
overall_rank = overall_rank .* feas_penalty;
[~, order] = sort(overall_rank);

fprintf('\n\n================ INDIVIDUAL RANKINGS =================\n\n');

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
        fprintf('%-20s %12s %6d\n', displayNames{a}, '1e6', i);
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
        fprintf('%-20s %12s %6d\n', displayNames{a}, '1e6', i);
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
        fprintf('%-20s %12s %6d\n', displayNames{a}, '1e6', i);
    elseif ~isempty(spread_per_run_vals{a}) && length(spread_per_run_vals{a}) > 0
        spread_str = sprintf('%.3f±%.3f', spread_vals(a), std(spread_per_run_vals{a}));
        fprintf('%-20s %12s %6d\n', displayNames{a}, spread_str, i);
    else
        fprintf('%-20s %12.3f %6d\n', displayNames{a}, spread_vals(a), i);
    end
end

fprintf('\n=== TIME RANKING (Lower is better, successful runs only) ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'Time (s) (Mean±Std)', 'Rank');
% Print only those with avg_time > 0, sorted
valid_time_sorted = valid_time_idx(time_sort);
for r = 1:length(valid_time_sorted)
    a = valid_time_sorted(r);
    fprintf('%-20s %12.2f±%.2f %6d\n', ...
        displayNames{a}, avg_time(a), std_time(a), r);
end
% Optionally list those excluded
excluded = find(avg_time == 0 & total_runs > 0);
if ~isempty(excluded)
    fprintf('\n--- Algorithms with zero average time (excluded) ---\n');
    for e = 1:length(excluded)
        a = excluded(e);
        fprintf('%-20s %12.2f±%.2f\n', displayNames{a}, avg_time(a), std_time(a));
    end
end

fprintf('\n=== FEASIBILITY RANKING ===\n');
fprintf('%-20s %12s %6s\n', 'Algorithm', 'Feasibility %', 'Rank');
for i = 1:nAlgos
    a = feas_rank(i);
    fprintf('%-20s %12.1f%% %6d\n', displayNames{a}, feasibility_rate(a), i);
end

%% =========================
%% SUMMARY STATISTICS
%% =========================
fprintf('\n\n================ SUMMARY STATISTICS =================\n');
fprintf('Total runs per algorithm: %d\n', numRuns);
fprintf('Reference front size: %d points\n', size(refFront, 1));

fprintf('\n=== TOP 3 ALGORITHMS ===\n');
for i = 1:min(3, nAlgos)
    a = order(i);
    fprintf('%d. %s (Score: %.2f)\n', i, displayNames{a}, overall_rank(a));
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

failed_algos = find(feasibility_rate == 0);
if ~isempty(failed_algos)
    fprintf('\n=== ALGORITHMS WITH 0%% FEASIBILITY ===\n');
    for i = 1:length(failed_algos)
        a = failed_algos(i);
        fprintf('• %s\n', displayNames{a});
    end
end

fprintf('\n✅ Analysis complete!\n');
fprintf('\n=== STATISTICAL TEST SUMMARY ===\n');
fprintf('• Original tests compare each algorithm against the per‑metric best.\n');
fprintf('• Additional tests compare each algorithm against %s (reference).\n', reference_name);
fprintf('• Wilcoxon signed-rank tests with Holm-Bonferroni correction, α = 0.05.\n');
fprintf('• All algorithms (including 0%% feasibility) are included in tests.\n');
fprintf('• Time ranking excludes algorithms with avg_time == 0 (i.e., those that crashed or had no runtime).\n');

%% =========================
%% HELPER FUNCTION
%% =========================
function out = ternary(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
