function build_PEM_Simulink()
%BUILD_PEM_SIMULINK  Create PEM_Plant.slx from scratch (robust, verbose).
%   Multi-strategy MATLAB Function script setter; forces compile between
%   stages so ports are correct before wiring. Idempotent.

    mdl = 'PEM_Plant';
    fprintf('\n[build_PEM_Simulink] Building %s.slx ...\n', mdl);

    if bdIsLoaded(mdl), close_system(mdl, 0); end
    if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
    if exist([mdl '.mdl'], 'file'), delete([mdl '.mdl']); end

    fprintf('  [1/6] Creating model skeleton...\n');
    new_system(mdl);
    load_system(mdl);

    set_param(mdl, ...
        'Solver',     'ode23t', ...
        'StopTime',   '200', ...
        'SaveOutput', 'on', ...
        'SaveFormat', 'Dataset', ...
        'SolverType', 'Variable-step');

    % ---------- Inputs ----------
    inTags  = {'SIM_A','SIM_N','SIM_J','SIM_T','SIM_P'};
    inNames = {'A_cell_cm2','N_cells','J_Acm2','T_C_in','P_bar'};
    ypos = [30 100 170 240 310];
    for k = 1:5
        add_block('simulink/Sources/Constant', [mdl '/' inNames{k}], ...
            'Value', inTags{k}, ...
            'Position', [30 ypos(k) 90 ypos(k)+30]);
    end

    fprintf('  [2/6] Adding MATLAB Function block...\n');
    add_block('simulink/User-Defined Functions/MATLAB Function', ...
        [mdl '/Electrochem'], 'Position', [200 60 380 300]);

    fprintf('  [3/6] Setting MATLAB Function script (multi-strategy)...\n');
    code = makeElectrochemCode();
    ok = setMATLABFunctionScript(mdl, 'Electrochem', code);

    if ~ok
        close_system(mdl, 0);
        error(['Could not set the MATLAB Function script on this MATLAB release.\n' ...
               '  Open %s.slx, double-click "Electrochem", paste the code from\n' ...
               '  makeElectrochemCode() in build_PEM_Simulink.m, then re-run.'], mdl);
    end

    % Force recompile so the block exposes its 5 in / 5 out ports
    fprintf('  [4/6] Compiling so port counts update...\n');
    try
        set_param(mdl, 'SimulationCommand', 'update');
    catch err
        fprintf('    (update warning, benign on first build: %s)\n', err.message);
    end

    % Verify ports
    numIn  = numel(find_system([mdl '/Electrochem'], 'FindAll','on', ...
                              'LookUnderMasks','all','BlockType','Inport'));
    numOut = numel(find_system([mdl '/Electrochem'], 'FindAll','on', ...
                              'LookUnderMasks','all','BlockType','Outport'));
    fprintf('    Electrochem ports: %d in / %d out\n', numIn, numOut);
    if numIn < 5 || numOut < 5
        close_system(mdl, 0);
        error(['MATLAB Function block still has only %d in / %d out ports.\n' ...
               '  Open %s.slx manually, open "Electrochem", and click Build.'], ...
               numIn, numOut, mdl);
    end

    fprintf('  [5/6] Wiring inputs, thermal loop, outputs...\n');

    % ---- 5 inputs -> Electrochem ----
    for k = 1:5
        add_line(mdl, sprintf('%s/1', inNames{k}), ...
                      sprintf('Electrochem/%d', k), 'autorouting', 'on');
    end

    % ---- Thermal ODE:  C dT/dt = P_heat - h_cool*(T - T_cool) ----
    add_block('simulink/Sources/Constant', [mdl '/T_cool_K'], ...
        'Value','333.15','Position',[430 400 490 430]);
    add_block('simulink/Sources/Constant', [mdl '/h_cool'], ...
        'Value','200','Position',[430 450 490 480]);
    add_block('simulink/Sources/Constant', [mdl '/C_th'], ...
        'Value','50000','Position',[430 500 490 530]);

    add_block('simulink/Signal Routing/From', [mdl '/From_T_K'], ...
        'GotoTag','T_K','Position',[430 340 480 360]);

    add_block('simulink/Math Operations/Sum', [mdl '/Sum_cool'], ...
        'Inputs','+-','Position',[520 345 550 375]);
    add_line(mdl,'From_T_K/1','Sum_cool/1','autorouting','on');
    add_line(mdl,'T_cool_K/1','Sum_cool/2','autorouting','on');

    add_block('simulink/Math Operations/Product', [mdl '/Prod_h'], ...
        'Inputs','2','Position',[590 350 630 380]);
    add_line(mdl,'Sum_cool/1','Prod_h/1','autorouting','on');
    add_line(mdl,'h_cool/1','Prod_h/2','autorouting','on');

    add_block('simulink/Math Operations/Sum', [mdl '/Sum_energy'], ...
        'Inputs','+-','Position',[670 200 700 240]);
    add_line(mdl,'Electrochem/5','Sum_energy/1','autorouting','on');
    add_line(mdl,'Prod_h/1','Sum_energy/2','autorouting','on');

    add_block('simulink/Math Operations/Math Function', [mdl '/Inv_Cth'], ...
        'Operator','reciprocal','Position',[630 500 670 530]);
    add_line(mdl,'C_th/1','Inv_Cth/1','autorouting','on');

    add_block('simulink/Math Operations/Product', [mdl '/Prod_invC'], ...
        'Inputs','2','Position',[740 205 780 235]);
    add_line(mdl,'Sum_energy/1','Prod_invC/1','autorouting','on');
    add_line(mdl,'Inv_Cth/1','Prod_invC/2','autorouting','on');

    add_block('simulink/Continuous/Integrator', [mdl '/Int_T'], ...
        'InitialCondition','333.15','Position',[820 200 860 240]);
    add_line(mdl,'Prod_invC/1','Int_T/1','autorouting','on');

    add_block('simulink/Signal Routing/Goto', [mdl '/Goto_T_K'], ...
        'GotoTag','T_K','Position',[900 210 940 230]);
    add_line(mdl,'Int_T/1','Goto_T_K/1','autorouting','on');

    % ---- Outputs ----
    add_block('simulink/Sinks/Out1',[mdl '/V_cell'],   'Position',[450  60 480  80]);
    add_line(mdl,'Electrochem/1','V_cell/1','autorouting','on');

    add_block('simulink/Sinks/Out1',[mdl '/P_stack'],  'Position',[450 120 480 140]);
    add_line(mdl,'Electrochem/3','P_stack/1','autorouting','on');

    add_block('simulink/Sinks/Out1',[mdl '/m_H2'],     'Position',[450 180 480 200]);
    add_line(mdl,'Electrochem/4','m_H2/1','autorouting','on');

    add_block('simulink/Sinks/Out1',[mdl '/T_cell_K'], 'Position',[980 210 1010 230]);
    add_line(mdl,'Int_T/1','T_cell_K/1','autorouting','on');

    fprintf('  [6/6] Saving...\n');
    save_system(mdl);
    fprintf('[build_PEM_Simulink] Done. Model: %s.slx\n\n', mdl);
end

% =====================================================================
function code = makeElectrochemCode()
    code = strjoin({ ...
        'function [V_cell, I_stack, P_stack, m_H2, P_heat] = fcn(A_cell_cm2, N_cells, J_Acm2, T_C, P_bar)'
        '    F      = 96485;'
        '    R      = 8.314;'
        '    M_H2   = 2.016e-3;'
        '    E_rev  = 1.42;'
        '    alpha  = 0.448938;'
        '    j0     = 5.428292e-02;'
        '    R_ohm  = 0.169475;'
        '    k_conc = 4e-6;'
        '    n_e    = 2;'
        '    T_K    = T_C + 273.15;'
        '    E_rev_T = E_rev - 0.0009*(T_C - 75);'
        '    f = R*T_K/(alpha*n_e*F);'
        '    eta_act  = f * asinh(J_Acm2/(2*j0));'
        '    eta_ohm  = R_ohm * J_Acm2;'
        '    eta_conc = k_conc * J_Acm2^2;'
        '    V_nernst = (R*T_K/(2*F)) * log(P_bar/1.0);'
        '    V_cell   = E_rev_T + eta_act + eta_ohm + eta_conc + V_nernst;'
        '    I_stack  = J_Acm2 * A_cell_cm2;'
        '    P_stack  = N_cells * V_cell * I_stack;'
        '    m_H2     = N_cells * I_stack / (2*F) * M_H2;'
        '    P_heat   = 0.25 * P_stack;'
        'end'}, newline);
end

% =====================================================================
function ok = setMATLABFunctionScript(mdl, blockPath, scriptText)
    fullPath = [mdl '/' blockPath];
    ok = false;

    % --- Strategy 0: force a compile so the inner chart exists ---
    try
        set_param(mdl, 'SimulationCommand', 'update');
        fprintf('    [setMFScript] pre-compile ok\n');
    catch err
        fprintf('    [setMFScript] pre-compile warning: %s\n', err.message);
    end

    % --- Strategy 1: get_param(...MATLABFunctionConfiguration).Script = ... ---
    try
        cfg = get_param(fullPath, 'MATLABFunctionConfiguration');
        cfg.Script = scriptText;
        fprintf('    [setMFScript] Strategy 1 (MATLABFunctionConfiguration) OK\n');
        ok = true;
        return;
    catch err1
        fprintf('    [setMFScript] Strategy 1 failed: %s\n', err1.message);
    end

    % --- Strategy 2: sfroot.find by path ---
    try
        rt = sfroot;
        charts = rt.find('-isa', 'Stateflow.EMChart');
        fprintf('    [setMFScript] Strategy 2: found %d EMChart(s)\n', numel(charts));
        for i = 1:numel(charts)
            p = charts(i).Path;
            fprintf('      chart(%d).Path = %s\n', i, p);
            % Simulink stores internally with backslash in some versions
            pNorm = strrep(p, '\', '/');
            if strcmpi(pNorm, fullPath) || endsWith(pNorm, ['/' blockPath])
                charts(i).Script = scriptText;
                fprintf('    [setMFScript] Strategy 2 OK on chart %d\n', i);
                ok = true;
                return;
            end
        end
        fprintf('    [setMFScript] Strategy 2: no matching chart path.\n');
    catch err2
        fprintf('    [setMFScript] Strategy 2 failed: %s\n', err2.message);
    end

    % --- Strategy 3: set_param on the config object ---
    try
        cfg = get_param(fullPath, 'MATLABFunctionConfiguration');
        set_param(cfg, 'Script', scriptText);
        fprintf('    [setMFScript] Strategy 3 (set_param on cfg) OK\n');
        ok = true;
        return;
    catch err3
        fprintf('    [setMFScript] Strategy 3 failed: %s\n', err3.message);
    end

    % --- Strategy 4: close/reopen to force chart creation, then retry S2 ---
    try
        fprintf('    [setMFScript] Strategy 4: close/reopen to force chart creation\n');
        save_system(mdl);
        close_system(mdl, 0);
        load_system(mdl);
        set_param(mdl, 'SimulationCommand', 'update');
        rt = sfroot;
        charts = rt.find('-isa', 'Stateflow.EMChart');
        for i = 1:numel(charts)
            pNorm = strrep(charts(i).Path, '\', '/');
            if endsWith(pNorm, ['/' blockPath])
                charts(i).Script = scriptText;
                fprintf('    [setMFScript] Strategy 4 OK on chart %d\n', i);
                ok = true;
                return;
            end
        end
    catch err4
        fprintf('    [setMFScript] Strategy 4 failed: %s\n', err4.message);
    end
end