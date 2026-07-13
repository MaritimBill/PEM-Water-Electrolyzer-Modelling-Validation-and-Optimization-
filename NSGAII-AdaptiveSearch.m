classdef NSGAII_AdaptiveSearch < ALGORITHM
% <multi> <real/integer> <constrained>
% NSGA-II with adaptive operators and local search (no change to selection).

    methods
        function main(Algorithm, Problem)
            maxGen = Problem.maxFE / Problem.N;
            pm_init = 1 / Problem.D;
            pm_final = 0.5 / Problem.D;
            eta_c_init = 15;
            eta_c_final = 30;
            localFreq = 5;  % local search every 5 generations

            Population = Problem.Initialization();
            [Population, FrontNo, CrowdDis] = EnvironmentalSelection(Population, Problem.N);

            while Algorithm.NotTerminated(Population)
                gen = floor((Problem.FE - Problem.N) / Problem.N) + 1;
                progress = gen / maxGen;

                % Adaptive parameters
                pm = pm_init + (pm_final - pm_init) * progress;
                eta_c = eta_c_init + (eta_c_final - eta_c_init) * progress;

                % Mating and offspring (standard NSGA‑II)
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                Offspring = OperatorGA(Problem, Population(MatingPool), ...
                                       {1, eta_c, pm, 20});

                % Environmental selection (ORIGINAL – no changes)
                [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);

                % --- Local search and diversity injection ---
                if mod(gen, localFreq) == 0
                    feasible = all(Population.cons <= 1e-6, 2);
                    if any(feasible)
                        % Best feasible by sum of objectives
                        [~, bestIdx] = min(sum(Population(feasible).objs, 2));
                        bestSol = Population(feasible(bestIdx));

                        % 3 random perturbations (2% of range)
                        for step = 1:3
                            perturb = 0.02 * (Problem.upper - Problem.lower) .* (2*rand(1, Problem.D) - 1);
                            newDec = bestSol.decs + perturb;
                            newDec = max(Problem.lower, min(Problem.upper, newDec));
                            newSol = Problem.Evaluation(newDec);
                            if all(newSol.cons <= 1e-6) && sum(newSol.objs) < sum(bestSol.objs)
                                bestSol = newSol;
                            end
                        end
                        % Replace worst (by crowding) in the last front
                        worstFront = find(FrontNo == max(FrontNo));
                        [~, worstRank] = sort(CrowdDis(worstFront));
                        Population(worstFront(worstRank(1))) = bestSol;
                    end

                    % Inject random feasible if diversity is low
                    validCrowd = CrowdDis(~isinf(CrowdDis) & ~isnan(CrowdDis));
                    if ~isempty(validCrowd) && std(validCrowd) < 0.05 * mean(validCrowd)
                        worstFront = find(FrontNo == max(FrontNo));
                        for k = 1:2
                            randDec = Problem.lower + rand(1, Problem.D) .* (Problem.upper - Problem.lower);
                            randSol = Problem.Evaluation(randDec);
                            if all(randSol.cons <= 1e-6)
                                [~, worstRank] = sort(CrowdDis(worstFront));
                                Population(worstFront(worstRank(1))) = randSol;
                            end
                        end
                    end
                end
            end
        end
    end
end
