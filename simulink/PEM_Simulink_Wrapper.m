classdef PEM_Simulink_Wrapper < handle
%PEM_SIMULINK_WRAPPER  Fast-restart wrapper for PEM_Plant.slx.
%
%   Usage:
%       w = PEM_Simulink_Wrapper();
%       [V, P, mH2, T_K] = w.simulate(A_cm2, N, J_Am2, T_C, P_Pa);
%
%   First call includes model compile (~seconds); subsequent calls
%   typically 5-20 ms thanks to Accelerator + FastRestart.

    properties (Access = private)
        mdl       = 'PEM_Plant';
        isLoaded  = false;
    end

    methods
        function obj = PEM_Simulink_Wrapper()
            if exist([obj.mdl '.slx'], 'file') == 0
                fprintf('PEM_Plant.slx not found -> building it now...\n');
                build_PEM_Simulink();
            end
            if ~bdIsLoaded(obj.mdl)
                load_system(obj.mdl);
            end
            set_param(obj.mdl, ...
                'SimulationMode', 'accelerator', ...
                'FastRestart',    'on');
            obj.isLoaded = true;
        end

        function [V_cell, P_stack, m_H2, T_K] = simulate(obj, A_cm2, N, J_Am2, T_C, P_Pa)
            if ~obj.isLoaded
                error('Wrapper not initialised.');
            end
            J_Acm2 = J_Am2 / 1e4;
            P_bar  = P_Pa  / 1e5;

            in = Simulink.SimulationInput(obj.mdl);
            in = in.setVariable('SIM_A', A_cm2);
            in = in.setVariable('SIM_N', N);
            in = in.setVariable('SIM_J', J_Acm2);
            in = in.setVariable('SIM_T', T_C);
            in = in.setVariable('SIM_P', P_bar);

            out = sim(in);

            V_cell  = out.yout{1}.Values.Data(end);
            P_stack = out.yout{2}.Values.Data(end);
            m_H2    = out.yout{3}.Values.Data(end);
            T_K     = out.yout{4}.Values.Data(end);
        end

        function delete(obj)
            try
                set_param(obj.mdl, 'FastRestart', 'off');
            catch
            end
            obj.isLoaded = false;
        end
    end
end