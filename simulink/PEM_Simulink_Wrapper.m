classdef PEM_Simulink_Wrapper < handle
%PEM_SIMULINK_WRAPPER  Fast wrapper around PEM_Plant.slx.
%
%   Construction:
%       w = PEM_Simulink_Wrapper('PEM_Plant');
%       % loads the model, sets accelerator mode, warms up the target
%
%   Simulation:
%       [V, P, m] = w.simulate(A, N, J, T_C, P_bar);
%
%   Inputs — PHYSICAL engineering units, NOT scaled decision variables:
%       A      = active area           [cm^2]    e.g.  50 ... 200
%       N      = number of cells       [-]       e.g.  10 ... 100
%       J      = current density       [A/cm^2]  e.g. 0.1 ... 2
%       T_C    = cell temperature      [degC]    e.g.  40 ... 90
%       P_bar  = H2 partial pressure   [bar]     e.g.   1 ... 30
%
%   Outputs — scalars from the last sample of the model run:
%       V      = single-cell voltage   [V]
%       P      = stack power           [W]
%       m      = H2 mass flow          [kg/s]
%
%   Notes
%   -----
%   * The model is left loaded between calls. Delete the wrapper
%     (clear w) if you need to reload PEM_Plant.slx.
%   * out.yout is a Dataset with three UNNAMED elements because the
%     model uses SaveOutput='on' (legacy output ports). Order matches
%     the Outport creation order in build_PEM_Simulink.m:
%         1 = V_cell,  2 = P_stack,  3 = m_H2

    properties (SetAccess = immutable)
        mdl = 'PEM_Plant';   % model name (no .slx extension)
    end

    methods
        % -----------------------------------------------------------
        function obj = PEM_Simulink_Wrapper(mdl)
            if nargin >= 1 && ~isempty(mdl)
                obj.mdl = mdl;
            end

            % ---- Load the model if it isn't already in memory ----
            if ~bdIsLoaded(obj.mdl)
                if exist([obj.mdl '.slx'], 'file') ~= 2 && ...
                   exist([obj.mdl '.mdl'], 'file') ~= 2
                    error('PEM_Simulink_Wrapper:ModelNotFound', ...
                          'Cannot find %s.slx on the MATLAB path. ' , ...
                          'Run build_PEM_Simulink first.', obj.mdl);
                end
                load_system(obj.mdl);
            end

            % ---- Configure run mode ----
            set_param(obj.mdl, ...
                'SimulationMode', 'accelerator', ...
                'StopTime',       '1');

            % ---- Warm-up compile so accelerator target exists ----
            obj.warm();
        end

        % -----------------------------------------------------------
        function warm(obj)
        %WARM  Compile the model once so the accelerator target is ready.
        %      Safe to call repeatedly; ignores the expected first-build
        %      warning about unconnected ports on a half-built model.
            try
                set_param(obj.mdl, 'SimulationCommand', 'update');
            catch
                % benign — ports are fine after build_PEM_Simulink finished
            end
        end

        % -----------------------------------------------------------
        function [V_cell, P_stack, m_H2] = simulate(obj, A, N, J, T_C, P_bar)
        %SIMULATE  One PEM_Plant evaluation; returns the three scalars.

            in = Simulink.SimulationInput(obj.mdl);
            in = in.setModelParameter('SimulationMode', 'accelerator', ...
                                      'StopTime',       '1');
            in = in.setVariable('SIM_A', A);
            in = in.setVariable('SIM_N', N);
            in = in.setVariable('SIM_J', J);
            in = in.setVariable('SIM_T', T_C);
            in = in.setVariable('SIM_P', P_bar);

            out = sim(in);

            % ---- Read the three output-port signals ----
            % out.yout: 3 unnamed Dataset elements, order fixed by the
            % Outport block creation order in build_PEM_Simulink.m.
            if out.yout.numElements < 3
                error('PEM_Simulink_Wrapper:NotEnoughOutputs', ...
                      'Model produced %d outputs, expected 3.', ...
                      out.yout.numElements);
            end

            V_cell  = out.yout{1}.Values.Data(end);
            P_stack = out.yout{2}.Values.Data(end);
            m_H2    = out.yout{3}.Values.Data(end);
        end

        % -----------------------------------------------------------
        function delete(obj)
        %DELETE  Called on clear / destruction.
        %        The model is intentionally left loaded for speed.
        %        Uncomment the line below if you want automatic unload.
            % if bdIsLoaded(obj.mdl), close_system(obj.mdl, 0); end
        end
    end
end
