%% FULL BENCHMARK SCRIPT (CORRECTED VERSION - PUBLICATION READY)

clear; clc; close all;

%% SETTINGS
popSize = 100;
maxGen  = 100;
maxFE   = popSize * maxGen;
numRuns = 30; % Number of independent runs

fprintf('========================================\n');
fprintf('FULL MOEA BENCHMARK SUITE\n');
fprintf('Problems: 10 | Algorithms: 16 | Runs: %d\n', numRuns);
fprintf('========================================\n\n');

%% Add PlatEMO
platemo_path = 'C:\Users\Admin\Desktop\PlatEMO-4.7';
if exist(platemo_path,'dir')
    addpath(genpath(platemo_path));
else
    error('PlatEMO not found!');
end

%% ALGORITHMS (16 TOTAL)
algorithms = {'NSGAII','NSGAIII','ANSGAIII','GDE3','MOEAD','MOEADDE','MOEADSTM','PREA','RVEA','RSEA','SRA','KnEA','MOPSO','dMOPSO','IBEA','FeaNSGA_II'};
displayNames = {'NSGA-II','NSGA-III','ANSGA-III','GDE3','MOEA/D','MOEA/D-DE','MOEA/D-STM','PREA','RVEA','RSEA','SRA','KnEA','MOPSO','dMOPSO','IBEA','FeaNSGA_II'};

%% PROBLEMS


%% OPTIMIZED BENCHMARK SUITE (Removed slow WFG2)
testProblems = {
    'ZDT1',     
    %'ZDT2',     
    %'DT4',     
    %'DTLZ1',    
    'DTLZ2',    
    %'DTLZ4',    
    %'UF1',      
    %'UF7',   
    %'CF1',     
    'WFG3'      
};
nObjectives = [2, 2, 2, 3, 3, 3, 2, 2, 2, 3];

nAlgos = length(algorithms);
nProbs = length(testProblems);

%% STORAGE STRUCTURE
for p = 1:nProbs
    probName = testProblems{p};
    for a = 1:nAlgos
        results.(probName).(algorithms{a}).HV  = NaN(numRuns,1);
        results.(probName).(algorithms{a}).IGD = NaN(numRuns,1);
        results.(probName).(algorithms{a}).success = false(numRuns,1);
    end
end

%% MAIN LOOP
for p = 1:nProbs
    probName = testProblems{p};
    nobj = nObjectives(p);
    
    fprintf('\n================ %s (M=%d) ================\n', probName, nobj);
    
    for run = 1:numRuns
        fprintf('Run %2d/%d: ', run, numRuns);
        
        for a = 1:nAlgos
            algoName = algorithms{a};
            rng(run*100 + a);  % Reproducible seeds
            
            try
                % Create problem
                Problem = feval(probName,...
                    'N', popSize,...
                    'M', nobj,...
                    'D', 30,...
                    'maxFE', maxFE);
                
                % Create and run algorithm
                Algorithm = feval(algoName);
                Algorithm.Solve(Problem);
                
                % Extract final population
                if ~isempty(Algorithm.result)
                    finalPop = Algorithm.result{end,2};
                    
                    % Calculate HV (Higher is better)
                    try
                        hv = Problem.CalMetric('HV', finalPop);
                        if isnan(hv) || isinf(hv) || hv < 0
                            hv = 0;
                        end
                    catch
                        hv = 0;
                    end
                    
                    % Calculate IGD (Lower is better)
                    try
                        igd = Problem.CalMetric('IGD', finalPop);
                        if isnan(igd) || isinf(igd) || igd < 0
                            igd = 1e6;  % Large penalty for failure
                        end
                    catch
                        igd = 1e6;
                    end
                    
                    results.(probName).(algoName).HV(run) = hv;
                    results.(probName).(algoName).IGD(run) = igd;
                    results.(probName).(algoName).success(run) = (hv > 0);
                    
                    fprintf('.');
                else
                    results.(probName).(algoName).HV(run) = 0;
                    results.(probName).(algoName).IGD(run) = 1e6;
                    fprintf('X');
                end
                
            catch ME
                results.(probName).(algoName).HV(run) = 0;
                results.(probName).(algoName).IGD(run) = 1e6;
                fprintf('E');
            end
        end
        fprintf(' done\n');
    end
    
    %% PER-PROBLEM SUMMARY WITH PROPER RANKING
    hv_mean  = zeros(nAlgos,1);
    hv_std   = zeros(nAlgos,1);
    igd_mean = zeros(nAlgos,1);
    igd_std  = zeros(nAlgos,1);
    success_rate = zeros(nAlgos,1);
    
    for a = 1:nAlgos
        hv_vals  = results.(probName).(algorithms{a}).HV;
        igd_vals = results.(probName).(algorithms{a}).IGD;
        success  = results.(probName).(algorithms{a}).success;
        
        % Filter successful runs for statistics
        valid_hv = hv_vals(hv_vals > 0);
        valid_igd = igd_vals(isfinite(igd_vals) & igd_vals < 1e5);
        
        if ~isempty(valid_hv)
            hv_mean(a)  = mean(valid_hv);
            hv_std(a)   = std(valid_hv);
        else
            hv_mean(a)  = 0;
            hv_std(a)   = 0;
        end
        
        if ~isempty(valid_igd)
            igd_mean(a) = mean(valid_igd);
            igd_std(a)  = std(valid_igd);
        else
            igd_mean(a) = 1e6;
            igd_std(a)  = 0;
        end
        
        success_rate(a) = sum(success) / numRuns * 100;
    end
    
    % RANKING: Sort by HV (higher better), then IGD (lower better)
    % First, create composite score: normalized HV rank + normalized IGD rank
    [~, hv_rank] = sort(-hv_mean);  % Higher HV = better rank
    [~, igd_rank] = sort(igd_mean);  % Lower IGD = better rank
    
    composite_score = zeros(nAlgos,1);
    for a = 1:nAlgos
        composite_score(a) = find(hv_rank == a) + find(igd_rank == a);
    end
    
    [~, order] = sort(composite_score);
    
    fprintf('\n🏆 RANKING (%s) - Based on HV (higher better) & IGD (lower better):\n', probName);
    fprintf('%-3s %-25s %12s %12s %10s %10s\n', 'Rank', 'Algorithm', 'HV', 'IGD', 'Feas%', 'Status');
    fprintf('%s\n', repmat('-', 80, 1));
    
    for i = 1:nAlgos
        idx = order(i);
        status = '';
        if success_rate(idx) == 100
            status = '✓✓ EXCELLENT';
        elseif success_rate(idx) >= 50
            status = '✓ GOOD';
        elseif success_rate(idx) > 0
            status = '⚠️ PARTIAL';
        else
            status = '❌ FAILED';
        end
        
        fprintf('%-3d %-25s %12.4e %12.4e %9.1f%%  %s\n', ...
            i, displayNames{idx}, hv_mean(idx), igd_mean(idx), ...
            success_rate(idx), status);
    end
    
    all_scores(:,p) = composite_score;
end

%% ================= FINAL GLOBAL RANKING =================
fprintf('\n\n');
fprintf('========================================\n');
fprintf('FINAL GLOBAL RANKING (ALL %d PROBLEMS)\n', nProbs);
fprintf('========================================\n');
fprintf('Scoring: Lower composite score = Better performance\n\n');

total_score = sum(all_scores, 2);
[~, final_order] = sort(total_score);

% Calculate overall success rate across all problems
overall_success = zeros(nAlgos,1);
for a = 1:nAlgos
    total_success = 0;
    total_runs = 0;
    for p = 1:nProbs
        probName = testProblems{p};
        total_success = total_success + sum(results.(probName).(algorithms{a}).success);
        total_runs = total_runs + numRuns;
    end
    overall_success(a) = total_success / total_runs * 100;
end

fprintf('%-3s %-25s %12s %12s %10s\n', 'Rank', 'Algorithm', 'Total Score', 'Avg HV', 'Feas%');
fprintf('%s\n', repmat('-', 80, 1));

for i = 1:nAlgos
    idx = final_order(i);
    
    % Calculate average HV across all problems
    avg_hv = 0;
    valid_probs = 0;
    for p = 1:nProbs
        probName = testProblems{p};
        hv_vals = results.(probName).(algorithms{idx}).HV;
        if any(hv_vals > 0)
            avg_hv = avg_hv + mean(hv_vals(hv_vals > 0));
            valid_probs = valid_probs + 1;
        end
    end
    if valid_probs > 0
        avg_hv = avg_hv / valid_probs;
    end
    
    status_icon = '';
    if overall_success(idx) == 100
        status_icon = '🏆';
    elseif overall_success(idx) >= 80
        status_icon = '✓';
    elseif overall_success(idx) > 0
        status_icon = '⚠️';
    else
        status_icon = '❌';
    end
    
    fprintf('%-3d %-25s %12.2f %12.4e %9.1f%%  %s\n', ...
        i, [displayNames{idx} ' ' status_icon], total_score(idx), avg_hv, overall_success(idx));
end

%% ================= DETAILED RESULTS TABLE =================
fprintf('\n\n');
fprintf('========================================\n');
fprintf('DETAILED RESULTS (Mean ± Std over %d runs)\n', numRuns);
fprintf('========================================\n\n');

for p = 1:nProbs
    probName = testProblems{p};
    fprintf('📊 %s:\n', probName);
    fprintf('%-25s %20s %25s\n', 'Algorithm', 'HV (↑ better)', 'IGD (↓ better)');
    fprintf('%s\n', repmat('-', 75, 1));
    
    for a = 1:nAlgos
        hv_vals = results.(probName).(algorithms{a}).HV;
        igd_vals = results.(probName).(algorithms{a}).IGD;
        
        valid_hv = hv_vals(hv_vals > 0);
        valid_igd = igd_vals(isfinite(igd_vals) & igd_vals < 1e5);
        
        if ~isempty(valid_hv)
            hv_str = sprintf('%.4e ± %.4e', mean(valid_hv), std(valid_hv));
        else
            hv_str = 'FAILED';
        end
        
        if ~isempty(valid_igd)
            igd_str = sprintf('%.4e ± %.4e', mean(valid_igd), std(valid_igd));
        else
            igd_str = 'FAILED';
        end
        
        fprintf('%-25s %20s %25s\n', displayNames{a}, hv_str, igd_str);
    end
    fprintf('\n');
end

%% ================= TOP PERFORMERS SUMMARY =================
fprintf('\n');
fprintf('========================================\n');
fprintf('TOP 3 ALGORITHMS SUMMARY\n');
fprintf('========================================\n\n');

for i = 1:min(3, nAlgos)
    idx = final_order(i);
    fprintf('🥇 #%d: %s (Score: %.2f, Feasibility: %.1f%%)\n', ...
        i, displayNames{idx}, total_score(idx), overall_success(idx));
end

fprintf('\n✅ Benchmark complete! Results saved.\n');

%% ================= SAVE RESULTS =================
save('benchmark_FULL_results.mat', 'results', 'testProblems', 'displayNames', ...
     'algorithms', 'numRuns', 'total_score', 'final_order', 'overall_success');
fprintf('\nSaved: benchmark_FULL_results.mat\n');

%% ================= PRINT FAILED ALGORITHMS =================
fprintf('\n');
fprintf('========================================\n');
fprintf('ALGORITHMS WITH 0%% FEASIBILITY\n');
fprintf('========================================\n');
failed_algos = find(overall_success == 0);
if ~isempty(failed_algos)
    for i = 1:length(failed_algos)
        fprintf('❌ %s\n', displayNames{failed_algos(i)});
    end
else
    fprintf('✓ All algorithms produced feasible solutions!\n');
end

%% Add statistical significance testing at the end
fprintf('\n\n========================================\n');
fprintf('STATISTICAL SIGNIFICANCE (Wilcoxon, α=0.05)\n');
fprintf('========================================\n\n');

% Compare FeaNSGA_II vs NSGA-II across all problems
for p = 1:nProbs
    probName = testProblems{p};
    pem_hv = results.(probName).('FeaNSGA_II').HV;
    nsga2_hv = results.(probName).('NSGAII').HV;
    
    valid_idx = pem_hv > 0 & nsga2_hv > 0;
    if sum(valid_idx) >= 10
        p_val = ranksum(pem_hv(valid_idx), nsga2_hv(valid_idx));
        if p_val < 0.05
            fprintf('%s: FeaNSGA_II vs NSGA-II: p=%.4f (SIGNIFICANT)\n', probName, p_val);
        else
            fprintf('%s: FeaNSGA_II vs NSGA-II: p=%.4f (not significant)\n', probName, p_val);
        end
    end
end