function build_PEM_Simulink()
%BUILD_PEM_SIMULINK  Create PEM_Plant.slx (electrochemistry only, 3 outputs).
%   Model I/O:
%       Inputs : A_cell_cm2, N_cells, J_Acm2, T_C, P_bar
%       Outputs: V_cell, P_stack, m_H2

    mdl = 'PEM_Plant';
    fprintf('\n[build_PEM_Simulink] Building %s.slx ...\n', mdl);

    if bdIsLoaded(mdl), close_system(mdl, 0); end
    if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
    if exist([mdl '.mdl'], 'file'), delete([mdl '.mdl']); end

    fprintf('  [1/5] Creating model skeleton...\n');
    new_system(mdl);
    load_system(mdl);

    set_param(mdl, ...
        'Solver',     'ode23t', ...
        'StopTime',   '1', ...
        'SaveOutput', 'on', ...
        'SaveFormat', 'Dataset', ...
        'SolverType', 'Variable-step');

    % ---------- Inputs (5) ----------
    inTags  = {'SIM_A','SIM_N','SIM_J','SIM_T','SIM_P'};
    inNames = {'A_cell_cm2','N_cells','J_Acm2','T_C_in','P_bar'};
    ypos = [30 100 170 240 310];
    for k = 1:5
        add_block('simulink/Sources/Constant', [mdl '/' inNames{k}], ...
            'Value', inTags{k}, ...
            'Position', [30 ypos(k) 90 ypos(k)+30]);
    end

    % ---------- MATLAB Function block ----------
    fprintf('  [2/5] Adding MATLAB Function block...\n');
    add_block('simulink/User-Defined Functions/MATLAB Function', ...
        [mdl '/Electrochem'], 'Position', [200 60 380 320]);

    fprintf('  [3/5] Setting MATLAB Function script...\n');
    code = makeElectrochemCode();
    ok = setMATLABFunctionScript(mdl, 'Electrochem', code);

    if ~ok
        close_system(mdl, 0);
        error(['Could not set the MATLAB Function script on this MATLAB release.\n' ...
               '  Open %s.slx, double-click "Electrochem", paste the code from\n' ...
               '  makeElectrochemCode() in build_PEM_Simulink.m, then re-run.'], mdl);
    end

    % Force recompile so the block exposes its 5 in / 3 out ports
    fprintf('  [4/5] Compiling so port counts update...\n');
    try
        set_param(mdl, 'SimulationCommand', 'update');
    catch err
        fprintf('    (update warning, benign on first build: %s)\n', err.message);
    end

    numIn  = numel(find_system([mdl '/Electrochem'], 'FindAll','on', ...
                              'LookUnderMasks','all','BlockType','Inport'));
    numOut = numel(find_system([mdl '/Electrochem'], 'FindAll','on', ...
                              'LookUnderMasks','all','BlockType','Outport'));
    fprintf('    Electrochem ports: %d in / %d out\n', numIn, numOut);
    if numIn < 5 || numOut < 3
        close_system(mdl, 0);
        error(['MATLAB Function block has %d in / %d out ports (need 5 / 3).\n' ...
               '  Open %s.slx, open "Electrochem", click Build.'], ...
               numIn, numOut, mdl);
    end

    % ---------- Wire inputs -> Electrochem ----------
    fprintf('  [5/5] Wiring inputs and outputs...\n');
    for k = 1:5
        add_line(mdl, sprintf('%s/1', inNames{k}), ...
                      sprintf('Electrochem/%d', k), 'autorouting', 'on');
    end

    % ---------- Outputs (3) ----------
    add_block('simulink/Sinks/Out1',[mdl '/V_cell'],  'Position',[450  80 480 100]);
    add_line(mdl,'Electrochem/1','V_cell/1','autorouting','on');

    add_block('simulink/Sinks/Out1',[mdl '/P_stack'], 'Position',[450 160 480 180]);
    add_line(mdl,'Electrochem/2','P_stack/1','autorouting','on');

    add_block('simulink/Sinks/Out1',[mdl '/m_H2'],    'Position',[450 240 480 260]);
    add_line(mdl,'Electrochem/3','m_H2/1','autorouting','on');

    save_system(mdl);
    fprintf('[build_PEM_Simulink] Done. Model: %s.slx\n\n', mdl);
end

% =====================================================================
function code = makeElectrochemCode()
    code = strjoin({ ...
        'function [V_cell, P_stack, m_H2] = fcn(A_cell_cm2, N_cells, J_Acm2, T_C, P_bar)'
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
        'end'}, newline);
end

% =====================================================================
function ok = setMATLABFunctionScript(mdl, blockPath, scriptText)
    fullPath = [mdl '/' blockPath];
    ok = false;

    % --- Strategy 0: force compile so the inner chart exists ---
    try
        set_param(mdl, 'SimulationCommand', 'update');
        fprintf('    [setMFScript] pre-compile ok\n');
    catch err
        fprintf('    [setMFScript] pre-compile warning: %s\n', err.message);
    end

    % --- Strategy 1: MATLABFunctionConfiguration ---
    try
        cfg = get_param(fullPath, 'MATLABFunctionConfiguration');
        cfg.Script = scriptText;
        fprintf('    [setMFScript] Strategy 1 OK\n');
        ok = true;
        return;
    catch err1
        fprintf('    [setMFScript] Strategy 1 failed: %s\n', err1.message);
    end

    % --- Strategy 2: sfroot.find EMChart by path ---
    try
        rt = sfroot;
        charts = rt.find('-isa', 'Stateflow.EMChart');
        for i = 1:numel(charts)
            pNorm = strrep(charts(i).Path, '\', '/');
            if strcmpi(pNorm, fullPath) || endsWith(pNorm, ['/' blockPath])
                charts(i).Script = scriptText;
                fprintf('    [setMFScript] Strategy 2 OK on chart %d\n', i);
                ok = true;
                return;
            end
        end
    catch err2
        fprintf('    [setMFScript] Strategy 2 failed: %s\n', err2.message);
    end

    % --- Strategy 3: set_param on config object ---
    try
        cfg = get_param(fullPath, 'MATLABFunctionConfiguration');
        set_param(cfg, 'Script', scriptText);
        fprintf('    [setMFScript] Strategy 3 OK\n');
        ok = true;
        return;
    catch err3
        fprintf('    [setMFScript] Strategy 3 failed: %s\n', err3.message);
    end

    % --- Strategy 4: close/reopen then retry S2 ---
    try
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
