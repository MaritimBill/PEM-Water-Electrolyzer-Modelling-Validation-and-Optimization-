%% ========================================================================
%% COMPREHENSIVE 19-ALGORITHM COMPARISON + SIMULINK DIGITAL-TWIN VALIDATION
%% File: comparisonss.m   (lives at PlatEMO root, next to Algorithms/)
%% ========================================================================
clear; clc; close all;

% ===================== PLATEMO PATH BOOTSTRAP =====================
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

platemoRoot = scriptDir;

if exist(fullfile(platemoRoot, 'Algorithms'), 'dir') ~= 7
    error(['Could not locate PlatEMO root. Expected:\n' ...
           '  %s\n' ...
           '  with an "Algorithms" subfolder inside it.\n' ...
           'Currently running from: %s'], platemoRoot, scriptDir);
end

addpath(genpath(platemoRoot));
fprintf('PlatEMO root on path: %s\n', platemoRoot);

assert(exist('ALGORITHM', 'class') == 8, ...
    'ALGORITHM class not on path.');
assert(exist('NSGAII', 'class') == 8, ...
    'NSGAII not found.');
assert(exist('PEM_Problem_Scaled', 'class') == 8, ...
    'PEM_Problem_Scaled not found.');
if exist('NSGAII_AdaptiveSearch', 'class') ~= 8
    warning(['NSGAII_AdaptiveSearch not found on path. ' ...
             'Is it in PlatEMO/Sims/?']);
end
if exist('PEM_Simulink_Wrapper', 'class') ~= 8
    warning(['PEM_Simulink_Wrapper not found on path. ' ...
             'Is it in PlatEMO/Sims/?']);
end
fprintf('Path OK.\n');
% ==================================================================

% ------------------------- CONFIG ----------------------------------------
numRuns    = 1;       % <<<<< SET TO 30 FOR FULL STATISTICAL ANALYSIS
popSize    = 100;
maxGen     = 100;
maxFE      = popSize * maxGen;
ptsPerAlgo = 20;      % <<<<< Simulink sims per algorithm (validation stage)
% -------------------------------------------------------------------------

algorithms = {'NSGAII','NSGAIII','ANSGAIII','GDE3','MOEAD','MOEADDE','MOEADSTM',...
              'ARMOEA','RVEA','RSEA','SRA','KnEA','BCEIBEA','SPEA2','IBEA',...
              'NSGAII_AdaptiveSearch','dMOPSO','MOPSO','PREA'};
displayNames = {'NSGA-II','NSGA-III','ANSGA-III','GDE3','MOEA/D','MOEA/D-DE',...
                'MOEA/D-STM','AR-MOEA','RVEA','RSEA','SRA','KnEA','BCE-IBEA',...
                'SPEA2','IBEA','NSGAII_AdaptiveSearch','dMOPSO','MOPSO','PREA'};

nAlgos = length(algorithms);

fprintf('\n================================================================\n');
fprintf('  COMPARISON START\n');
fprintf('  Algorithms : %d\n', nAlgos);
fprintf('  Runs/algo  : %d\n', numRuns);
fprintf('  Pop / Gen  : %d / %d  (FE/run = %d)\n', popSize, maxGen, maxFE);
if numRuns < 5
    fprintf('  NOTE: numRuns<5 -> Wilcoxon tests will be SKIPPED.\n');
end
fprintf('================================================================\n');

% Storage
hv_per_run     = cell(nAlgos, 1);
igd_per_run    = cell(nAlgos, 1);
gd_per_run     = cell(nAlgos, 1);
spread_per_run = cell(nAlgos, 1);
time_per_run   = cell(nAlgos, 1);
feas_per_run   = cell(nAlgos, 1);

results = struct;
all_feasible_objs = [];

feasible_count   = zeros(nAlgos, 1);
total_runs       = zeros(nAlgos, 1);
total_time       = zeros(nAlgos, 1);
feasibility_rate = zeros(nAlgos, 1);

%% ========================================================================
%% RUN
%% ========================================================================
tAllStart = tic;
for run = 1:numRuns
    fprintf('\n========== RUN %d / %d ==========\n', run, numRuns);
    for a = 1:nAlgos
        algoName = algorithms{a};
        rng(run * 100 + a, 'twister');

        Problem = PEM_Problem_Scaled();
        Problem.N = popSize;
        Problem.maxFE = maxFE;

        pct = 100 * ((run-1)*nAlgos + (a-1)) / (numRuns*nAlgos);
        elapsed = toc(tAllStart);
        if pct > 0
            etaSec = elapsed * (100 - pct) / pct;
            etaStr = sprintf('%d:%02d', floor(etaSec/60), mod(round(etaSec),60));
        else
            etaStr = '--:--';
        end
        fprintf('[%5.1f%% | ETA %s] %-20s ', pct, etaStr, displayNames{a});

        if ~isfield(results, algoName)
            results.(algoName).objs           = [];
            results.(algoName).decs           = [];
            results.(algoName).time           = [];
            results.(algoName).feasible_count = 0;
            results.(algoName).per_run_objs   = {};
            results.(algoName).per_run_decs   = {};
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
                fprintf('X (no population, %.2fs)\n', runtime);
                results.(algoName).time(end + 1)         = runtime;
                results.(algoName).per_run_objs{end + 1} = [];
                results.(algoName).per_run_decs{end + 1} = [];
                total_runs(a) = total_runs(a) + 1;
                total_time(a) = total_time(a) + runtime;
                continue;
            end

            objs     = cat(1, finalPop.objs);
            cons     = cat(1, finalPop.cons);
            decs     = cat(1, finalPop.decs);
            feasible = all(cons <= 0, 2);
            n_feasible = sum(feasible);

            results.(algoName).time(end + 1) = runtime;

            if n_feasible > 0
                objs_f = objs(feasible, :);
                decs_f = decs(feasible, :);

                all_feasible_objs = [all_feasible_objs; objs_f];              %#ok<AGROW>
                results.(algoName).objs = [results.(algoName).objs; objs_f];  %#ok<AGROW>
                results.(algoName).decs = [results.(algoName).decs; decs_f];  %#ok<AGROW>
                results.(algoName).feasible_count = ...
                    results.(algoName).feasible_count + n_feasible;
                results.(algoName).per_run_objs{end + 1} = objs_f;
                results.(algoName).per_run_decs{end + 1} = decs_f;

                feasible_count(a) = feasible_count(a) + n_feasible;
                fprintf('ok (%d feasible, %.2fs)\n', n_feasible, runtime);
            else
                results.(algoName).per_run_objs{end + 1} = [];
                results.(algoName).per_run_decs{end + 1} = [];
                fprintf('no feasible (%.2fs)\n', runtime);
            end

            total_runs(a) = total_runs(a) + 1;
            total_time(a) = total_time(a) + runtime;

        catch err
            fprintf('ERROR: %s\n', err.message);
            if ~isempty(err.stack)
                fprintf('         at %s (line %d)\n', ...
                    err.stack(1).name, err.stack(1).line);
            end
            if ~isfield(results, algoName)
                results.(algoName).objs           = [];
                results.(algoName).decs           = [];
                results.(algoName).time           = [];
                results.(algoName).feasible_count = 0;
                results.(algoName).per_run_objs   = {};
                results.(algoName).per_run_decs   = {};
            end
            results.(algoName).time(end + 1)         = 0;
            results.(algoName).per_run_objs{end + 1} = [];
            results.(algoName).per_run_decs{end + 1} = [];
            total_runs(a) = total_runs(a) + 1;
        end
    end

    fprintf('\n----- RUN %d SUMMARY -----\n', run);
    for a2 = 1:nAlgos
        name = algorithms{a2};
        if isfield(results, name) && ~isempty(results.(name).per_run_objs)
            last = results.(name).per_run_objs{end};
            if ~isempty(last)
                fprintf('   %-20s  %d feasible\n', displayNames{a2}, size(last,1));
            end
        end
    end
end

for a = 1:nAlgos
    if total_runs(a) > 0
        expected_total = total_runs(a) * 100;
        feasibility_rate(a) = feasible_count(a) / expected_total * 100;
    else
        feasibility_rate(a) = 0;
    end
end

%% ========================================================================
%% REFERENCE FRONT
%% ========================================================================
if isempty(all_feasible_objs)
    error('No feasible solutions found across any algorithm.');
end

ND = NDSort(all_feasible_objs, size(all_feasible_objs, 1)) == 1;
refFront = all_feasible_objs(ND, :);

minF = min(refFront, [], 1);
maxF = max(refFront, [], 1);
refFront_norm = (refFront - minF) ./ (maxF - minF + 1e-12);

fprintf('\n[Reference front] %d Pareto points.\n', size(refFront,1));

%% ========================================================================
%% METRICS
%% ========================================================================
hv_vals     = inf(nAlgos, 1);
igd_vals    = inf(nAlgos, 1);
gd_vals     = inf(nAlgos, 1);
spread_vals = inf(nAlgos, 1);
avg_time    = zeros(nAlgos, 1);
std_time    = zeros(nAlgos, 1);

hv_per_run_vals     = cell(nAlgos, 1);
igd_per_run_vals    = cell(nAlgos, 1);
gd_per_run_vals     = cell(nAlgos, 1);
spread_per_run_vals = cell(nAlgos, 1);
time_per_run_vals   = cell(nAlgos, 1);
feas_per_run_vals   = cell(nAlgos, 1);

refPoint = ones(1, 3) * 1.1;

for a = 1:nAlgos
    algoName = algorithms{a};

    if ~isfield(results, algoName) || isempty(results.(algoName).objs)
        hv_vals(a) = 0; igd_vals(a) = 1e6; gd_vals(a) = 1e6; spread_vals(a) = 1e6;
        avg_time(a) = total_time(a) / max(1, total_runs(a));
        std_time(a) = 0;
        hv_per_run_vals{a} = []; igd_per_run_vals{a} = []; gd_per_run_vals{a} = [];
        spread_per_run_vals{a} = []; time_per_run_vals{a} = []; feas_per_run_vals{a} = [];
        continue;
    end

    times = results.(algoName).time;
    avg_time(a) = mean(times);
    std_time(a) = std(times);

    feas_per_run = zeros(total_runs(a), 1);
    for r = 1:total_runs(a)
        run_objs = results.(algoName).per_run_objs{r};
        if ~isempty(run_objs)
            feas_per_run(r) = size(run_objs, 1);
        end
    end
    feas_per_run_vals{a} = feas_per_run;

    hv_run     = zeros(total_runs(a), 1);
    igd_run    = zeros(total_runs(a), 1);
    gd_run     = zeros(total_runs(a), 1);
    spread_run = zeros(total_runs(a), 1);

    for r = 1:total_runs(a)
        run_objs = results.(algoName).per_run_objs{r};

        if isempty(run_objs)
            hv_run(r) = 0; igd_run(r) = 1e6; gd_run(r) = 1e6; spread_run(r) = 1e6;
            continue;
        end

        objs_n = (run_objs - minF) ./ (maxF - minF + 1e-12);

        Nsample = 10000;
        samples = rand(Nsample, 3) .* refPoint;
        dominated = false(Nsample, 1);
        for i = 1:size(objs_n, 1)
            dominated = dominated | all(samples >= objs_n(i, :), 2);
        end
        hv_run(r) = sum(dominated) / Nsample * prod(refPoint);

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

    hv_per_run_vals{a}     = hv_run;
    igd_per_run_vals{a}    = igd_run;
    gd_per_run_vals{a}     = gd_run;
    spread_per_run_vals{a} = spread_run;
    time_per_run_vals{a}   = times;

    hv_vals(a)     = mean(hv_run);
    igd_vals(a)    = mean(igd_run);
    gd_vals(a)     = mean(gd_run);
    spread_vals(a) = mean(spread_run);
end

%% ========================================================================
%% STATISTICAL TESTS
%% ========================================================================
fprintf('\n\n================ STATISTICAL ANALYSIS =================\n');
metrics_names   = {'HV', 'IGD', 'GD', 'Spread', 'Time'};
metrics_vals    = {hv_vals, igd_vals, gd_vals, spread_vals, avg_time};
metrics_per_run = {hv_per_run_vals, igd_per_run_vals, gd_per_run_vals, ...
                   spread_per_run_vals, time_per_run_vals};
better_direction = {'higher', 'lower', 'lower', 'lower', 'lower'};

if numRuns < 5
    fprintf('Skipped: numRuns = %d (< 5).\n', numRuns);
    fprintf('Increase numRuns at the top of comparisonss.m for Wilcoxon tests.\n');
else
    for m = 1:length(metrics_names)
        metric_name    = metrics_names{m};
        metric_vals    = metrics_vals{m};
        metric_per_run = metrics_per_run{m};
        better         = better_direction{m};

        fprintf('\n=== %s Metric (Wilcoxon signed-rank, Holm-Bonferroni) ===\n', metric_name);
        valid_idx = find(~isinf(metric_vals) & ~isnan(metric_vals));
        if isempty(valid_idx)
            fprintf('No valid metric values.\n'); continue;
        end

        if strcmp(better, 'higher')
            [best_val, best_pos] = max(metric_vals(valid_idx));
        else
            [best_val, best_pos] = min(metric_vals(valid_idx));
        end
        best_metric_idx = valid_idx(best_pos);
        fprintf('Best: %s (%.4f)\n\n', displayNames{best_metric_idx}, best_val);

        p_values = []; algo_indices = [];
        for a = 1:nAlgos
            if a == best_metric_idx, continue; end
            if isempty(metric_per_run{a}) || isempty(metric_per_run{best_metric_idx})
                continue;
            end
            n_runs = min(length(metric_per_run{a}), length(metric_per_run{best_metric_idx}));
            if n_runs < 3, continue; end
            data_best = metric_per_run{best_metric_idx}(1:n_runs);
            data_a    = metric_per_run{a}(1:n_runs);
            data_best = data_best(~isnan(data_best) & ~isinf(data_best));
            data_a    = data_a(~isnan(data_a)    & ~isinf(data_a));
            if length(data_best) < 3 || length(data_a) < 3, continue; end
            n = min(length(data_best), length(data_a));
            data_best = data_best(1:n); data_a = data_a(1:n);
            if all(data_a == data_best), continue; end
            try
                if strcmp(better, 'higher')
                    p = signrank(data_a, data_best);
                else
                    p = signrank(data_best, data_a);
                end
                p_values = [p_values; p];        %#ok<AGROW>
                algo_indices = [algo_indices; a]; %#ok<AGROW>
            catch
            end
        end

        if ~isempty(p_values)
            [sorted_p, sort_idx] = sort(p_values);
            n_tests = length(sorted_p);
            alpha = 0.05;
            fprintf('%-20s %12s %12s %s\n', 'Algorithm', 'p-value', 'Adj.Thresh', 'Sig');
            fprintf('%s\n', repmat('-', 60, 1));
            for k = 1:n_tests
                adj_threshold = alpha / (n_tests - k + 1);
                is_sig = sorted_p(k) < adj_threshold;
                a_idx  = algo_indices(sort_idx(k));
                if is_sig, sigStr = 'YES'; else, sigStr = 'NO'; end
                fprintf('%-20s %12.4e %12.4f %s\n', ...
                    displayNames{a_idx}, sorted_p(k), adj_threshold, sigStr);
            end
        else
            fprintf('No valid comparisons.\n');
        end
    end

    reference_name = 'NSGAII_AdaptiveSearch';
    ref_idx = find(strcmp(algorithms, reference_name));
    if ~isempty(ref_idx)
        fprintf('\n\n================ vs %s =================\n', reference_name);
        for m = 1:length(metrics_names)
            metric_name    = metrics_names{m};
            metric_per_run = metrics_per_run{m};
            better         = better_direction{m};
            fprintf('\n=== %s ===\n', metric_name);
            fprintf('%-20s %12s %12s %s\n', 'Algorithm', 'p-value', 'Adj.Thresh', 'Sig');
            fprintf('%s\n', repmat('-', 60, 1));

            data_ref = metric_per_run{ref_idx};
            data_ref = data_ref(~isnan(data_ref) & ~isinf(data_ref));
            if isempty(data_ref)
                fprintf('Reference has no valid data.\n'); continue;
            end

            p_values = []; algo_indices = [];
            for a = 1:nAlgos
                if a == ref_idx, continue; end
                data_a = metric_per_run{a};
                data_a = data_a(~isnan(data_a) & ~isinf(data_a));
                if isempty(data_a), continue; end
                n = min(length(data_ref), length(data_a));
                if n < 3, continue; end
                x = data_ref(1:n); y = data_a(1:n);
                if all(x == y), continue; end
                try
                    if strcmp(better, 'higher')
                        p = signrank(y, x);
                    else
                        p = signrank(x, y);
                    end
                    p_values = [p_values; p];         %#ok<AGROW>
                    algo_indices = [algo_indices; a]; %#ok<AGROW>
                catch
                end
            end

            if ~isempty(p_values)
                [sorted_p, sort_idx] = sort(p_values);
                n_tests = length(sorted_p);
                alpha = 0.05;
                for k = 1:n_tests
                    adj_threshold = alpha / (n_tests - k + 1);
                    is_sig = sorted_p(k) < adj_threshold;
                    a_idx  = algo_indices(sort_idx(k));
                    if is_sig, sigStr = 'YES'; else, sigStr = 'NO'; end
                    fprintf('%-20s %12.4e %12.4f %s\n', ...
                        displayNames{a_idx}, sorted_p(k), adj_threshold, sigStr);
                end
            else
                fprintf('No valid comparisons.\n');
            end
        end
    end
end

%% ========================================================================
%% SIMULINK DIGITAL-TWIN VALIDATION (per-algorithm, live progress)
%%
%%  *** CORRECTED SECTION ***
%%  The wrapper expects PHYSICAL engineering units [A cm^2, N cells,
%%  J A/cm^2, T degC, P bar].  The decision variables in `d` are SCALED
%%  to [0,1]^5 by PEM_Problem_Scaled.  We un-scale them here using
%%  unscaleDecision() before calling SLW.simulate().
%% ========================================================================
fprintf('\n\n================ SIMULINK DIGITAL-TWIN VALIDATION =================\n');

SLW = [];
try
    if exist('PEM_Plant.slx', 'file') ~= 4
        build_PEM_Simulink();
    end
    SLW = PEM_Simulink_Wrapper();
    fprintf('Simulink model loaded (accelerator + fast restart).\n');
catch err
    fprintf('Simulink unavailable: %s\n', err.message);
    fprintf('Skipping digital-twin validation.\n');
end

if ~isempty(SLW)
    P = PEM_Problem_Scaled();

    algo_valid_meanV = nan(nAlgos,1);
    algo_valid_meanP = nan(nAlgos,1);
    algo_valid_meanM = nan(nAlgos,1);
    algo_valid_maxP  = nan(nAlgos,1);
    algo_valid_nPts  = zeros(nAlgos,1);
    algo_valid_time  = zeros(nAlgos,1);

    % Print the un-scaling mapping once so the log shows what units
    % are being fed to the wrapper.
    [~, ~, info] = unscaleDecision(P, zeros(1,5));   %#ok<ASGLU>
    fprintf('  Un-scaling:  phys = lb + d .* (ub - lb)\n');
    fprintf('     lb = [%s]\n', num2str(info.lb, '%g '));
    fprintf('     ub = [%s]\n', num2str(info.ub, '%g '));
    fprintf('  (order: A_cm2, N_cells, J_Acm2, T_degC, P_bar)\n\n');

    tAllSL = tic;
    for a = 1:nAlgos
        algoName = algorithms{a};

        if ~isfield(results, algoName) || isempty(results.(algoName).objs)
            fprintf('\n>> %-22s  skipped (0 feasible)\n', displayNames{a});
            continue;
        end

        algo_objs = results.(algoName).objs;
        algo_decs = results.(algoName).decs;

        ND_a   = NDSort(algo_objs, size(algo_objs,1)) == 1;
        pf_dec = algo_decs(ND_a, :);
        nPf    = size(pf_dec, 1);

        if nPf == 0
            fprintf('\n>> %-22s  skipped (empty Pareto front)\n', displayNames{a});
            continue;
        end

        nChk = min(nPf, ptsPerAlgo);
        if nPf > ptsPerAlgo
            chkIdx = round(linspace(1, nPf, ptsPerAlgo));
        else
            chkIdx = 1:nPf;
        end

        header = sprintf('\n>> %-22s  %d feasible | %d Pareto | validating %d', ...
            displayNames{a}, size(algo_objs,1), nPf, nChk);
        fprintf('%s\n', header);

        relV = nan(nChk,1); relP = nan(nChk,1); relM = nan(nChk,1);

        tA = tic;
        for k = 1:nChk
            d = pf_dec(chkIdx(k), :);           % scaled, [0,1]^5
            m = P.CalculateMetrics(d);          % internal un-scaling

            % ---- Un-scale d to physical engineering units ----
            d_phys = unscaleDecision(P, d);

            try
                [V_sl, P_sl, m_sl] = SLW.simulate( ...
                    d_phys(1), d_phys(2), d_phys(3), d_phys(4), d_phys(5));

                relV(k) = 100 * abs(V_sl - m.V_cell  ) / max(abs(m.V_cell),   1);
                relP(k) = 100 * abs(P_sl - m.P_stack ) / max(abs(m.P_stack),  1);
                relM(k) = 100 * abs(m_sl - m.mass_H2 ) / max(abs(m.mass_H2),  1e-12);

                line = sprintf(['   [%2d/%2d]  phys=[%7.3f %6.2f %6.3f %6.2f %6.2f]' ...
                                '  V=%.4f V  P=%9.2f W  mH2=%.3e kg/s' ...
                                '  |  dV=%.4f%%  dP=%.4f%%  dm=%.4f%%'], ...
                    k, nChk, d_phys(1), d_phys(2), d_phys(3), d_phys(4), d_phys(5), ...
                    V_sl, P_sl, m_sl, relV(k), relP(k), relM(k));
                fprintf('%s\n', line);
            catch err
                line = sprintf('   [%2d/%2d]  SL error: %s', k, nChk, err.message);
                fprintf('%s\n', line);
            end
        end
        dt_a = toc(tA);

        algo_valid_meanV(a) = mean(relV, 'omitnan');
        algo_valid_meanP(a) = mean(relP, 'omitnan');
        algo_valid_meanM(a) = mean(relM, 'omitnan');
        algo_valid_maxP(a)  = max(relP,  [], 'omitnan');
        algo_valid_nPts(a)  = nChk;
        algo_valid_time(a)  = dt_a;

        line = sprintf('   --- mean |dV|=%.4f%%  |dP|=%.4f%%  |dm|=%.4f%%  max |dP|=%.4f%%', ...
            algo_valid_meanV(a), algo_valid_meanP(a), algo_valid_meanM(a), algo_valid_maxP(a));
        fprintf('%s\n', line);

        line = sprintf('   --- wall %.2f s', dt_a);
        fprintf('%s\n', line);
    end
    dtAllSL = toc(tAllSL);

    fprintf('\n\n--- VALIDATION SUMMARY ---\n');
    fprintf('%-22s %5s %10s %10s %10s\n', ...
        'Algorithm', 'Pts', 'mean dV', 'mean dP', 'mean dm');
    fprintf('%s\n', repmat('-', 66, 1));
    for a = 1:nAlgos
        if algo_valid_nPts(a) > 0
            line = sprintf('%-22s %5d %10.4f %10.4f %10.4f', ...
                displayNames{a}, algo_valid_nPts(a), ...
                algo_valid_meanV(a), algo_valid_meanP(a), algo_valid_meanM(a));
            fprintf('%s\n', line);
        end
    end
    fprintf('\nTotal Simulink wall time: %.2f s\n', dtAllSL);

    delete(SLW);
end

%% ========================================================================
%% FINAL TABLE
%% ========================================================================
fprintf('\n\n================ FINAL =================\n');
fprintf('%-15s %16s %16s %16s %16s %14s %8s\n', ...
    'Algorithm', 'HV', 'IGD', 'GD', 'Spr', 'Time(s)', 'Feas%%');
fprintf('%s\n', repmat('-', 115, 1));

for i = 1:nAlgos
    a = i;
    feas_str = sprintf('%.1f%%', feasibility_rate(a));

    if ~isempty(hv_per_run_vals{a}) && numel(hv_per_run_vals{a}) > 1
        hv_str     = sprintf('%.4f+-%.4f', hv_vals(a),     std(hv_per_run_vals{a}));
        igd_str    = sprintf('%.4f+-%.4f', igd_vals(a),    std(igd_per_run_vals{a}));
        gd_str     = sprintf('%.4f+-%.4f', gd_vals(a),     std(gd_per_run_vals{a}));
        spread_str = sprintf('%.4f+-%.4f', spread_vals(a), std(spread_per_run_vals{a}));
        time_str   = sprintf('%.2f+-%.2f', avg_time(a),    std_time(a));
    else
        hv_str     = sprintf('%.4f', hv_vals(a));
        igd_str    = sprintf('%.4f', igd_vals(a));
        gd_str     = sprintf('%.4f', gd_vals(a));
        spread_str = sprintf('%.4f', spread_vals(a));
        time_str   = sprintf('%.2f', avg_time(a));
    end

    fprintf('%-15s %16s %16s %16s %16s %14s %8s\n', ...
        displayNames{a}, hv_str, igd_str, gd_str, spread_str, time_str, feas_str);
end

%% ========================================================================
%% RANKINGS
%% ========================================================================
[~, hv_rank]     = sort(hv_vals, 'descend');
[~, igd_rank]    = sort(igd_vals, 'ascend');
[~, gd_rank]     = sort(gd_vals, 'ascend');
[~, spread_rank] = sort(spread_vals, 'ascend');

valid_time_idx = find(avg_time > 0 & total_runs > 0);
time_sort = [];
if ~isempty(valid_time_idx)
    [~, time_sort] = sort(avg_time(valid_time_idx));
    time_rank = inf(nAlgos, 1);
    for r = 1:length(valid_time_idx)
        time_rank(valid_time_idx(time_sort(r))) = r;
    end
else
    time_rank = inf(nAlgos, 1);
end
time_rank(isinf(time_rank)) = nAlgos + 1;

[~, feas_rank] = sort(feasibility_rate, 'descend');

rankHV = zeros(nAlgos,1); rankIGD = zeros(nAlgos,1); rankGD = zeros(nAlgos,1);
rankSpread = zeros(nAlgos,1); rankFeas = zeros(nAlgos,1);
for i = 1:nAlgos
    rankHV(hv_rank(i))         = i;
    rankIGD(igd_rank(i))       = i;
    rankGD(gd_rank(i))         = i;
    rankSpread(spread_rank(i)) = i;
    rankFeas(feas_rank(i))     = i;
end

feas_penalty = ones(nAlgos, 1);
for a = 1:nAlgos
    if feasibility_rate(a) < 50
        feas_penalty(a) = 1e6;
    elseif feasibility_rate(a) < 90
        feas_penalty(a) = 10;
    end
end

overall_rank = 0.35*rankHV + 0.20*rankIGD + 0.10*rankGD + ...
               0.10*rankSpread + 0.10*time_rank + 0.15*rankFeas;
overall_rank = overall_rank .* feas_penalty;
[~, order] = sort(overall_rank);

fprintf('\n================ INDIVIDUAL RANKINGS =================\n');

fprintf('\n--- HV (higher better) ---\n');
for i = 1:nAlgos
    a = hv_rank(i);
    if ~isempty(hv_per_run_vals{a}) && numel(hv_per_run_vals{a}) > 1
        s = sprintf('%.4f+-%.4f', hv_vals(a), std(hv_per_run_vals{a}));
    else
        s = sprintf('%.4f', hv_vals(a));
    end
    fprintf('%-20s %18s   #%d\n', displayNames{a}, s, i);
end

fprintf('\n--- IGD (lower better) ---\n');
for i = 1:nAlgos
    a = igd_rank(i);
    if igd_vals(a) >= 1e5
        s = '1e6';
    elseif ~isempty(igd_per_run_vals{a}) && numel(igd_per_run_vals{a}) > 1
        s = sprintf('%.4f+-%.4f', igd_vals(a), std(igd_per_run_vals{a}));
    else
        s = sprintf('%.4f', igd_vals(a));
    end
    fprintf('%-20s %18s   #%d\n', displayNames{a}, s, i);
end

fprintf('\n--- GD (lower better) ---\n');
for i = 1:nAlgos
    a = gd_rank(i);
    if gd_vals(a) >= 1e5
        s = '1e6';
    elseif ~isempty(gd_per_run_vals{a}) && numel(gd_per_run_vals{a}) > 1
        s = sprintf('%.4f+-%.4f', gd_vals(a), std(gd_per_run_vals{a}));
    else
        s = sprintf('%.4f', gd_vals(a));
    end
    fprintf('%-20s %18s   #%d\n', displayNames{a}, s, i);
end

fprintf('\n--- SPREAD (lower better) ---\n');
for i = 1:nAlgos
    a = spread_rank(i);
    if spread_vals(a) >= 1e5
        s = '1e6';
    elseif ~isempty(spread_per_run_vals{a}) && numel(spread_per_run_vals{a}) > 1
        s = sprintf('%.4f+-%.4f', spread_vals(a), std(spread_per_run_vals{a}));
    else
        s = sprintf('%.4f', spread_vals(a));
    end
    fprintf('%-20s %18s   #%d\n', displayNames{a}, s, i);
end

fprintf('\n--- TIME (lower better, ran-only) ---\n');
if ~isempty(time_sort)
    valid_time_sorted = valid_time_idx(time_sort);
    for r = 1:length(valid_time_sorted)
        a = valid_time_sorted(r);
        fprintf('%-20s %10.2f+-%.2f s   #%d\n', ...
            displayNames{a}, avg_time(a), std_time(a), r);
    end
else
    fprintf('  (no algorithms with avg_time > 0)\n');
end

fprintf('\n--- FEASIBILITY (higher better) ---\n');
for i = 1:nAlgos
    a = feas_rank(i);
    fprintf('%-20s %8.1f%%   #%d\n', displayNames{a}, feasibility_rate(a), i);
end

%% ========================================================================
%% SUMMARY
%% ========================================================================
fprintf('\n\n================ SUMMARY =================\n');
fprintf('Runs/algo          : %d\n', numRuns);
fprintf('Reference front    : %d points\n', size(refFront, 1));
fprintf('Total wall time    : %.2f s\n', toc(tAllStart));

fprintf('\n=== TOP 3 (overall rank) ===\n');
for i = 1:min(3, nAlgos)
    a = order(i);
    fprintf('%d. %-20s  score=%.2f\n', i, displayNames{a}, overall_rank(a));
    fprintf('     HV=%.4f  IGD=%.4f  GD=%.4f  Spr=%.4f  Time=%.2fs  Feas=%.1f%%\n', ...
        hv_vals(a), igd_vals(a), gd_vals(a), spread_vals(a), ...
        avg_time(a), feasibility_rate(a));
end

failed_algos = find(feasibility_rate == 0);
if ~isempty(failed_algos)
    fprintf('\n=== ZERO FEASIBILITY ===\n');
    for i = 1:length(failed_algos)
        fprintf('  - %s\n', displayNames{failed_algos(i)});
    end
end

fprintf('\nDone.\n');

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================
function [d_phys, ok, info] = unscaleDecision(P, d)
%UNSCALEDECISION  Map a scaled decision vector in [0,1]^5 back to
%   physical engineering units [A_cm2, N_cells, J_Acm2, T_degC, P_bar].
%
%   The problem class PEM_Problem_Scaled stores its bounds in one of
%   several possible forms.  This helper tries them in order:
%       1. P.lower / P.upper           (PlatEMO convention)
%       2. P.lb    / P.ub
%       3. P.A_min/P.A_max, P.N_min/P.N_max, ... (per-variable fields)
%       4. [P.A_min ...] via a static method P.GetBounds()
%       5. Fallback: assume d is already physical (with a warning).
%
%   Outputs:
%       d_phys  physical decision vector (1x5)
%       ok      true if un-scaling was performed from known bounds
%       info    struct with fields .lb, .ub, .mode  (for logging)

    d = d(:)';
    info = struct('lb', [], 'ub', [], 'mode', '');

    % ---- Strategy 1: lower / upper ----
    if isprop(P, 'lower') && isprop(P, 'upper') && ...
       numel(P.lower) >= 5 && numel(P.upper) >= 5
        lb = P.lower(1:5);  ub = P.upper(1:5);
        d_phys = lb + d .* (ub - lb);
        info.lb = lb; info.ub = ub; info.mode = 'lower/upper';
        ok = true; return;
    end

    % ---- Strategy 2: lb / ub ----
    if isprop(P, 'lb') && isprop(P, 'ub') && ...
       numel(P.lb) >= 5 && numel(P.ub) >= 5
        lb = P.lb(1:5);  ub = P.ub(1:5);
        d_phys = lb + d .* (ub - lb);
        info.lb = lb; info.ub = ub; info.mode = 'lb/ub';
        ok = true; return;
    end

    % ---- Strategy 3: per-variable _min / _max ----
    names = {'A','N','J','T','P'};
    lb = zeros(1,5); ub = zeros(1,5); found = true;
    for i = 1:5
        nm = names{i};
        if isprop(P, [nm '_min']) && isprop(P, [nm '_max'])
            lb(i) = P.([nm '_min']);
            ub(i) = P.([nm '_max']);
        else
            found = false; break;
        end
    end
    if found
        d_phys = lb + d .* (ub - lb);
        info.lb = lb; info.ub = ub; info.mode = 'per-variable';
        ok = true; return;
    end

    % ---- Strategy 4: P.GetBounds() static method ----
    try
        [lb, ub] = P.GetBounds();
        if numel(lb) >= 5 && numel(ub) >= 5
            lb = lb(1:5); ub = ub(1:5);
            d_phys = lb + d .* (ub - lb);
            info.lb = lb; info.ub = ub; info.mode = 'GetBounds()';
            ok = true; return;
        end
    catch
    end

    % ---- Fallback: does d already look physical? ----
    % If any entry is outside [0,1] with meaningful margin, assume physical.
    if any(d < -1e-6) || any(d > 1 + 1e-6)
        warning('comparisonss:unscaleFallback', ...
                ['Could not find bounds on PEM_Problem_Scaled; ' ...
                 'assuming d is already physical.']);
        d_phys = d;
        info.mode = 'assumed-physical';
        ok = false;
        return;
    end

    % Nothing worked, and d looks scaled → cannot proceed
    error('comparisonss:unscaleFailed', ...
          ['Cannot determine how to un-scale the decision vector. ' ...
           'Add explicit un-scaling in comparisonss.m or expose bounds ' ...
           'on PEM_Problem_Scaled (lower/upper, lb/ub, or A_min/A_max...).']);
end
