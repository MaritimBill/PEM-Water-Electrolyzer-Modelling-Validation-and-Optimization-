classdef PEM_Problem_Scaled < PROBLEM
% PEM Water Electrolyzer Optimization Problem - Based on Current Electrolyzer Data
% Uses validated 4-parameter model from PEM_System_Calibrated_Current.m
%   x1: Cell area [cm²]              (10 - 500)
%   x2: Number of cells [-]          (5 - 100)
%   x3: Current density [A/m²]       (1000 - 20000) -> converts to 0.1-2.0 A/cm²
%   x4: Temperature [°C]             (50 - 85)
%   x5: Pressure [Pa]                (1e5 - 30e5)
%
% Calibration based on experimental data from PEMData.mat
% Model: V_cell = E_rev + (RT/(αnF))·asinh(J/(2J₀)) + R_ohm·J + k_conc·J²
% 4-parameter model with literature-derived bounds
% Combined RMSE = 0.018987 V, R² = 0.986711 (59 points)
%
% Objectives (MINIMIZE):
%   f1: Specific energy consumption [kWh/kg H2]
%   f2: Total active area [m²] with pressure scaling (cost proxy)
%   f3: Membrane stress index (durability proxy)
%
% Constraints:
%   - Cell voltage limits (1.3-2.5 V)
%   - Stack power limits (1-300 kW)
%   - Current density limits
%   - Temperature bounds
%   - Pressure bounds
%   - Minimum H2 production rate

    properties(Constant)
        %% ========== Physical Constants ==========
        F_const = 96485;    % Faraday constant [C/mol]
        R_gas   = 8.314;    % Gas constant [J/(mol·K)]
        M_H2    = 2.016e-3; % Molar mass of H2 [kg/mol]
        
        %% ========== Calibrated PEM Model Parameters ==========
        % 4-parameter model with literature-derived bounds
        % Calibrated using fmincon with k_conc ∈ [0, 0.02]
        % Combined RMSE = 0.018987 V, R² = 0.986711 (59 points)
        %
        % Literature-bounds test proved k_conc is negligible:
        %   - Unconstrained (allowed negative): k_conc = -0.0139 (unphysical)
        %   - Positive bounds [0, 0.1]: α=0.9, j₀=1e-3 (overfitting)
        %   - Extreme bounds [0, 1.0]: k_conc → 0
        %   - Literature bounds [0, 0.02]: k_conc = 4e-6 (effectively zero)
        %
        E_rev = 1.42;               % Reversible potential at 75°C [V]
        alpha = 0.448938;           % Transfer coefficient [-]
        j0 = 5.428292e-02;          % Exchange current density [A/cm²]
        R_ohm = 0.169475;           % Area-specific resistance [Ω·cm²]
        k_conc = 4e-6;              % Concentration coefficient [V/(A/cm²)²]
                                    % Retained as it converged from literature bounds
        n_e = 2;                    % Number of electrons per reaction
        
        %% ========== REALISTIC BOUNDS ==========
        LB = [10,    5,   1000,   50,     1e5];   % Lower bounds
        UB = [500, 100,  20000,   85,    30e5];   % Upper bounds
        
        %% ========== Design Requirements ==========
        %MIN_H2_RATE = 1.3577e-4;  % kg/s H2 (from optimization constraints)
        MIN_H2_RATE = 1.6971e-4;  % kg/s H2 
        MIN_POWER   = 1000;       % Minimum stack power [W]
        MAX_POWER   = 300e3;      % Maximum stack power [W]
        MAX_VOLTAGE = 2.5;        % Maximum cell voltage [V]
        MIN_VOLTAGE = 1.3;        % Minimum cell voltage [V]
        
        %% ========== Cost Model Coefficients ==========
        % cost = k1 * N * A + k2 * (A * P)^0.7 + k3 * N * P
        k1 = 150;      % Cell material cost scaling [$/cell·cm²]
        k2 = 0.1;      % Containment cost scaling [$/cm²·bar^0.7]
        k3 = 0.5;      % Pressure sealing cost scaling [$/cell·bar]
        
        %% ========== Reference Values for Normalization ==========
        REF_ENERGY = 50;        % Reference specific energy [kWh/kg]
        REF_AREA   = 1;         % Reference area [m²]
        REF_COST   = 15000;     % Reference cost [USD] - for normalization
        REF_STRESS = 1;         % Reference stress index
    end
    
    methods
        %% ========== Problem Definition ==========
        function Setting(obj)
            obj.M = 3;                          % Number of objectives
            obj.D = 5;                          % Number of decision variables
            obj.lower = obj.LB;                 % Lower bounds
            obj.upper = obj.UB;                 % Upper bounds
            obj.encoding = ones(1, obj.D);      % Real encoding
        end

        %% ========== Calculate Objective Values ==========
        function PopObj = CalObj(obj, PopDec)
            N = size(PopDec, 1);
            PopObj = zeros(N, obj.M);
            for i = 1:N
                F = obj.EvaluatePEM(PopDec(i, :));
                PopObj(i,:) = F;
            end
        end
        
        %% ========== Calculate Constraint Violations ==========
        function PopCon = CalCon(obj, PopDec)
            N = size(PopDec, 1);
            PopCon = zeros(N, 10);
            for i = 1:N
                [~, G] = obj.EvaluatePEM(PopDec(i, :));
                PopCon(i,:) = G;
            end
        end
        
        %% ========== Core PEM Evaluation (4-Parameter Model) ==========
        function [F, G] = EvaluatePEM(obj, x)
            %% STRICT BOUND ENFORCEMENT
            invalid_mask = isnan(x) | isinf(x) | (x < 0);
            x(invalid_mask) = obj.LB(invalid_mask);
            
            below_LB = x < obj.LB;
            above_UB = x > obj.UB;
            bound_violation = sum(abs(x - max(obj.LB, min(obj.UB, x))));
            
            if any(below_LB) || any(above_UB)
                F = [1e10, 1e10, 1e10];
                G = ones(1, 10) * 1e6;
                return;
            end
            
            x_clamped = max(obj.LB, min(obj.UB, x));
            
            %% Extract Decision Variables
            A_cell_cm2 = x_clamped(1);   % Cell area [cm²]
            N_cells = x_clamped(2);       % Number of cells [-]
            J_Am2 = x_clamped(3);        % Current density [A/m²]
            T_C = x_clamped(4);           % Temperature [°C]
            P_Pa = x_clamped(5);          % Pressure [Pa]
            
            %% Unit Conversions
            J_Acm2 = J_Am2 / 10000;       % Convert to [A/cm²]
            T_K = T_C + 273.15;           % Convert to [K]
            A_cell_m2 = A_cell_cm2 / 10000; % Convert to [m²]
            P_bar = P_Pa / 1e5;            % Convert to [bar]
            
            %% ========== 4-PARAMETER PEM MODEL ==========
            % V_cell = E_rev + (RT/(αnF))·asinh(J/(2J₀)) + R_ohm·J + k_conc·J²
            % Butler-Volmer kinetics with literature-bounded k_conc
            
            % Temperature-dependent reversible potential
            % E_rev decreases ~0.9 mV/°C
            E_rev_T = obj.E_rev - 0.0009 * (T_C - 75);
            
            % Activation overpotential (Butler-Volmer kinetics)
            % η_act = (RT/(αnF)) * asinh(J/(2j₀))
            if J_Acm2 > 0
                f = obj.R_gas * T_K / (obj.alpha * obj.n_e * obj.F_const);
                eta_act = f * asinh(J_Acm2 / (2 * obj.j0));
            else
                eta_act = 0;
            end
            
            % Ohmic overpotential
            eta_ohm = obj.R_ohm * J_Acm2;
            
            % Concentration overpotential (k_conc ≈ 4e-6)
            eta_conc = obj.k_conc * J_Acm2^2;
            
            % Total cell voltage
            V_cell = E_rev_T + eta_act + eta_ohm + eta_conc;
            
            % Pressure effect (Nernst equation)
            if P_Pa > 0
                P_ref = 1e5;  % 1 bar
                V_cell = V_cell + (obj.R_gas * T_K / (2 * obj.F_const)) * log(P_Pa / P_ref);
            end
            
            %% Production Calculations
            I = J_Am2 * A_cell_m2;              % Stack current [A]
            P_total = N_cells * V_cell * I;     % Stack power [W]
            
            % Hydrogen production (Faraday's law)
            n_H2 = (N_cells * I) / (2 * obj.F_const);  % [mol/s]
            mass_H2 = n_H2 * obj.M_H2;                  % [kg/s]
            
            %% ========== OBJECTIVES (MINIMIZE) ==========
            
            % Objective 1: Specific Energy Consumption [kWh/kg H2]
            if mass_H2 > 1e-12
                specific_energy = (P_total / 1000) / (mass_H2 * 3600);
                specific_energy = max(30, min(80, specific_energy));
            else
                specific_energy = 200;
            end
            F(1) = specific_energy;
            
            % Objective 2: Total Active Area with Pressure Scaling [m²] (Cost Proxy)
            A_cell_cm2 = A_cell_m2 * 10000;
            
            % Three cost components
            cell_material_cost = obj.k1 * N_cells * A_cell_cm2;
            containment_cost = obj.k2 * (A_cell_cm2 * P_bar) ^ 0.7;
            pressure_seal_cost = obj.k3 * N_cells * P_bar;
            total_cost = cell_material_cost + containment_cost + pressure_seal_cost;
            
            % Convert cost to [m²-equivalent]
            cost_per_m2 = 15000;
            total_area_equivalent = total_cost / cost_per_m2;
            F(2) = total_area_equivalent;
            
            % Objective 3: Enhanced Stress/Durability Proxy
            J_norm = J_Acm2 / 2.0;
            T_norm = (T_C - 50) / 35;
            
            T_ref_deg = 333.15;
            T_deg_factor = exp(5000 * (1/T_ref_deg - 1/T_K));
            J_deg_factor = J_norm^1.5;
            
            stress_index = (J_deg_factor + T_deg_factor) / 2;
            stress_index = max(0.1, min(2, stress_index));
            F(3) = stress_index;
            
            %% ========== CONSTRAINTS (violation > 0 = infeasible) ==========
            G = zeros(1, 10);
            
            % Voltage constraints
            G(1) = max(0, V_cell - obj.MAX_VOLTAGE);
            G(2) = max(0, obj.MIN_VOLTAGE - V_cell);
            
            % Power constraints
            G(3) = max(0, (P_total / 1000) - (obj.MAX_POWER / 1000));
            G(4) = max(0, obj.MIN_POWER - P_total);
            
            % Current density constraints
            G(5) = max(0, 500 - J_Am2);
            G(6) = max(0, J_Am2 - 20000);
            
            % Temperature constraints
            G(7) = max(0, T_C - 85) + max(0, 50 - T_C);
            
            % Pressure constraints
            G(8) = max(0, P_Pa - 30e5) + max(0, 1e5 - P_Pa);
            
            % Minimum H2 production constraint
            G(9) = max(0, obj.MIN_H2_RATE - mass_H2);
            
            % Bound violation penalty
            G(10) = bound_violation * 1000;
        end
        
        %% ========== Helper: Get Polarization Data for Plotting ==========
        function [J_range, V_range] = GetPolarizationCurve(obj, T_C)
            J_range = linspace(0.01, 2.0, 100)';
            T_K = T_C + 273.15;
            E_rev_T = obj.E_rev - 0.0009 * (T_C - 75);
            
            V_range = zeros(size(J_range));
            for i = 1:length(J_range)
                J = J_range(i);
                if J > 0
                    f = obj.R_gas * T_K / (obj.alpha * obj.n_e * obj.F_const);
                    eta_act = f * asinh(J / (2 * obj.j0));
                else
                    eta_act = 0;
                end
                eta_ohm = obj.R_ohm * J;
                eta_conc = obj.k_conc * J^2;
                V_range(i) = E_rev_T + eta_act + eta_ohm + eta_conc;
            end
        end
        
        %% ========== Helper: Calculate Metrics for a Design ==========
        function metrics = CalculateMetrics(obj, x)
            [F, G] = obj.EvaluatePEM(x);
            
            A_cell_cm2 = x(1);
            N_cells = x(2);
            J_Am2 = x(3);
            T_C = x(4);
            P_Pa = x(5);
            
            J_Acm2 = J_Am2 / 10000;
            A_cell_m2 = A_cell_cm2 / 10000;
            T_K = T_C + 273.15;
            
            E_rev_T = obj.E_rev - 0.0009 * (T_C - 75);
            if J_Acm2 > 0
                f = obj.R_gas * T_K / (obj.alpha * obj.n_e * obj.F_const);
                eta_act = f * asinh(J_Acm2 / (2 * obj.j0));
            else
                eta_act = 0;
            end
            eta_ohm = obj.R_ohm * J_Acm2;
            eta_conc = obj.k_conc * J_Acm2^2;
            V_cell = E_rev_T + eta_act + eta_ohm + eta_conc;
            
            if P_Pa > 0
                P_ref = 1e5;
                V_cell = V_cell + (obj.R_gas * T_K / (2 * obj.F_const)) * log(P_Pa / P_ref);
            end
            
            I = J_Am2 * A_cell_m2;
            P_stack = N_cells * V_cell * I;
            
            n_H2 = (N_cells * I) / (2 * obj.F_const);
            mass_H2 = n_H2 * obj.M_H2;
            
            metrics.V_cell = V_cell;
            metrics.I = I;
            metrics.P_stack = P_stack;
            metrics.mass_H2 = mass_H2;
            metrics.specific_energy = F(1);
            metrics.total_area = F(2);
            metrics.stress_index = F(3);
            metrics.constraints = G;
            metrics.is_feasible = all(G <= 1e-6);
        end
    end
end
