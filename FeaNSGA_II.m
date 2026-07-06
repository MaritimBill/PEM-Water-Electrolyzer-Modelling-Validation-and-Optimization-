classdef FeaNSGA_II < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% Feasibility-first NSGA-II variant with elite archive, periodic immigrants, and integer repair.

    methods
        function main(Algorithm,Problem)
            N = Problem.N;
            D = Problem.D;

            %% Initial population
            Population = Problem.Initialization();
            Population = FeaNSGA_II.RepairAndEvaluate(Problem,Population);
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N);

            %% Feasible nondominated archive
            Archive = FeaNSGA_II.GetFeasibleNondominated(Population);
            if isempty(Archive)
                Archive = Population;
            end

            gen = 0;

            while Algorithm.NotTerminated(Population)
                %% Constraint violation of current population
                CV = FeaNSGA_II.ConstraintViolation(Population);

                %% Feasibility-first tournament selection
                MatingPool = zeros(1,N);
                for i = 1:N
                    a = randi(N);
                    b = randi(N);
                    MatingPool(i) = FeaNSGA_II.FeasibilityTournament(a,b,FrontNo,CrowdDis,CV);
                end

                Parents = Population(MatingPool);

                %% Inject archive elites into the mating pool
                if ~isempty(Archive)
                    eliteNum  = min(max(2,round(0.15*N)), numel(Archive));
                    eliteIdx  = randi(numel(Archive),1,eliteNum);
                    replaceIdx = randperm(N,eliteNum);
                    Parents(replaceIdx) = Archive(eliteIdx);
                end

                %% Standard GA offspring
                Offspring = OperatorGA(Problem,Parents);
                Offspring = FeaNSGA_II.RepairAndEvaluate(Problem,Offspring);

                %% Periodic random immigrants to preserve diversity
                if mod(gen,10) == 0
                    k = max(2,round(0.08*N));
                    ImmDec = repmat(Problem.lower,k,1) + ...
                             rand(k,D).*repmat(Problem.upper-Problem.lower,k,1);

                    if D >= 2
                        ImmDec(:,2) = round(ImmDec(:,2));   % integer variable
                    end

                    ImmPop = Problem.Evaluation(ImmDec);
                    Offspring = [Offspring, ImmPop];
                end

                %% Environmental selection
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],N);

                %% Refresh archive
                NewArchive = FeaNSGA_II.GetFeasibleNondominated(Population);
                if ~isempty(NewArchive)
                    Archive = FeaNSGA_II.MergeArchive(Archive,NewArchive,N);
                end

                gen = gen + 1;
            end
        end
    end

    methods (Static, Access = private)

        function Pop = RepairAndEvaluate(Problem,Pop)
            if isempty(Pop)
                return;
            end

            Dec = Pop.decs;

            % Clip to bounds
            Dec = max(repmat(Problem.lower,size(Dec,1),1), ...
                      min(repmat(Problem.upper,size(Dec,1),1), Dec));

            % x2 = number of cells should be integer
            if size(Dec,2) >= 2
                Dec(:,2) = round(Dec(:,2));
                Dec(:,2) = max(Problem.lower(2), min(Problem.upper(2), Dec(:,2)));
            end

            % Repair invalid rows with random feasible samples
            bad = any(~isfinite(Dec),2);
            if any(bad)
                m = sum(bad);
                Dec(bad,:) = repmat(Problem.lower,m,1) + ...
                             rand(m,Problem.D).*repmat(Problem.upper-Problem.lower,m,1);
                if size(Dec,2) >= 2
                    Dec(bad,2) = round(Dec(bad,2));
                end
            end

            Pop = Problem.Evaluation(Dec);
        end

        function CV = ConstraintViolation(Pop)
            if isempty(Pop)
                CV = [];
                return;
            end

            Con = Pop.cons;
            if isempty(Con)
                CV = zeros(numel(Pop),1);
            else
                CV = sum(max(0,Con),2);
            end
        end

        function selected = FeasibilityTournament(idxA,idxB,FrontNo,CrowdDis,CV)
            feasA = CV(idxA) <= 1e-12;
            feasB = CV(idxB) <= 1e-12;

            % Feasible solutions always beat infeasible ones
            if feasA && ~feasB
                selected = idxA;
                return;
            elseif ~feasA && feasB
                selected = idxB;
                return;
            end

            % If both feasible: lower FrontNo first, then higher crowding
            if feasA && feasB
                if FrontNo(idxA) < FrontNo(idxB)
                    selected = idxA;
                elseif FrontNo(idxB) < FrontNo(idxA)
                    selected = idxB;
                else
                    if CrowdDis(idxA) > CrowdDis(idxB)
                        selected = idxA;
                    elseif CrowdDis(idxB) > CrowdDis(idxA)
                        selected = idxB;
                    else
                        if rand < 0.5
                            selected = idxA;
                        else
                            selected = idxB;
                        end
                    end
                end
            else
                % If both infeasible: lower constraint violation wins
                if CV(idxA) < CV(idxB)
                    selected = idxA;
                elseif CV(idxB) < CV(idxA)
                    selected = idxB;
                else
                    if rand < 0.5
                        selected = idxA;
                    else
                        selected = idxB;
                    end
                end
            end
        end

        function Feas = IsFeasible(Pop)
            if isempty(Pop)
                Feas = false(0,1);
                return;
            end

            Con = Pop.cons;
            if isempty(Con)
                Feas = true(numel(Pop),1);
            else
                Feas = all(Con <= 1e-12, 2);
            end
        end

        function Archive = GetFeasibleNondominated(Pop)
            Feas = FeaNSGA_II.IsFeasible(Pop);
            Archive = Pop(Feas);

            if isempty(Archive)
                return;
            end

            [FrontNo,~] = NDSort(Archive.objs,1);
            Archive = Archive(FrontNo == 1);
        end

        function Archive = MergeArchive(A,B,N)
            Archive = [A,B];
            if isempty(Archive)
                return;
            end

            Feas = FeaNSGA_II.IsFeasible(Archive);
            Archive = Archive(Feas);

            if isempty(Archive)
                return;
            end

            [FrontNo,~] = NDSort(Archive.objs,1);
            Archive = Archive(FrontNo == 1);

            % Keep archive bounded
            cap = max(3, 5*N);
            if numel(Archive) > cap
                idx = randperm(numel(Archive),cap);
                Archive = Archive(idx);
            end
        end
    end
end