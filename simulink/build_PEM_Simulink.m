function build_PEM_Simulink()
%BUILD_PEM_SIMULINK  Create PEM_Plant.slx (single MATLAB Function block).
%   Model I/O:
%     Inputs : A_cell_cm2, N_cells, J_Acm2, T_C_in, P_bar
%     Outputs: V_cell, P_stack, m_H2   (in this Outport order)

    mdl = 'PEM_Plant';
    fprintf('\n[build_PEM_Simulink] Building %s.slx ...\n', mdl);

    % ---- Defensive cleanup: Fast Restart can block close_system ----
    if bdIsLoaded(mdl)
        try, set_param(mdl, 'FastRestart', 'off');            catch, end
        try, set_param(mdl, 'SimulationCommand', 'stop');     catch, end
        try, set_param(mdl, 'SimulationCommand', 'terminate');catch, end
        try, close_system(mdl, 0);                            catch, end
    end
    if bdIsLoaded(mdl)
        try, bdclose(mdl); catch, end
    end
    if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
    if exist([mdl '.mdl'], 'file'), delete([mdl '.mdl']); end

    % ---- Seed base-workspace defaults BEFORE building the model ----
    % The Constant blocks reference SIM_A ... SIM_P.  Without these in
    % the base workspace, the block Parameters are "Invalid setting".
    if ~exist('SIM_A','var'), assignin('base','SIM_A', 100); end
    if ~exist('SIM_N','var'), assignin('base','SIM_N',  20); end
    if ~exist('SIM_J','var'), assignin('base','SIM_J',   1); end
    if ~exist('SIM_T','var'), assignin('base','SIM_T',  70); end
    if ~exist('SIM_P','var'), assignin('base','SIM_P',  10); end

    % ---- Skeleton ----
    fprintf('  [1/4] Creating model skeleton...\n');
    new_system(mdl);
    load_system(mdl);
    set_param(mdl, ...
        'Solver',     'ode23t', ...
        'StopTime',   '1', ...
        'SaveOutput', 'on', ...
        'SaveFormat', 'Dataset', ...
        'SolverType', 'Variable-step');

    % ---- Geometry (constants and outports centred on block ports) ----
    IN_X   =  40;  IN_W  = 80;  IN_H  = 30;
    EC_L   = 220;  EC_T  = 60;
    EC_W   = 200;  EC_H  = 300;
    OUT_X  = 500;  OUT_W = 30;  OUT_H = 14;

    inY  = EC_T + EC_H * (1:5) / 6;   % 110,160,210,260,310
    outY = EC_T + EC_H * (1:3) / 4;   % 135,210,285

    inTags  = {'SIM_A','SIM_N','SIM_J','SIM_T','SIM_P'};
    inNames = {'A_cell_cm2','N_cells','J_Acm2','T_C_in','P_bar'};
    outNames = {'V_cell','P_stack','m_H2'};

    % ---- Inputs (5) ----
    fprintf('  [2/4] Placing inputs...\n');
    for k = 1:5
        add_block('simulink/Sources/Constant', [mdl '/' inNames{k}], ...
            'Value',    inTags{k}, ...
            'Position', [IN_X, round(inY(k)-IN_H/2), ...
                         IN_X+IN_W, round(inY(k)+IN_H/2)]);
    end

    % ---- MATLAB Function block ----
    fprintf('  [3/4] Adding MATLAB Function block...\n');
    add_block('simulink/User-Defined Functions/MATLAB Function', ...
        [mdl '/Electrochem'], ...
        'Position', [EC_L, EC_T, EC_L+EC_W, EC_T+EC_H]);

    code = makeElectrochemCode();
    ok = setMATLABFunctionScript(mdl, 'Electrochem', code);
    if ~ok
        close_system(mdl, 0);
        error(['Could not set the MATLAB Function script.\n' ...
               'Open %s.slx, double-click "Electrochem", paste the code\n' ...
               'from makeElectrochemCode() in build_PEM_Simulink.m, then re-run.'], mdl);
    end

    % ---- Outports (built before wiring so add_line finds all ports) ----
    for j = 1:3
        add_block('simulink/Sinks/Out1', [mdl '/' outNames{j}], ...
            'Position', [OUT_X, round(outY(j)-OUT_H/2), ...
                         OUT_X+OUT_W, round(outY(j)+OUT_H/2)]);
    end

    % ---- Wire inputs + outputs (all straight horizontal segments) ----
    fprintf('  [4/4] Wiring inputs and outputs...\n');
    for k = 1:5
        add_line(mdl, sprintf('%s/1', inNames{k}), ...
                      sprintf('Electrochem/%d', k), 'autorouting','on');
    end
    for j = 1:3
        add_line(mdl, sprintf('Electrochem/%d', j), ...
                      sprintf('%s/1', outNames{j}), 'autorouting','on');
    end

    % ---- Now update: workspace is seeded and model is wired -> clean ----
    try
        set_param(mdl, 'SimulationCommand', 'update');
    catch err
        fprintf('    [build] update warning: %s\n', err.message);
        for i = 1:numel(err.cause)
            fprintf('      cause %d: %s\n', i, err.cause{i}.message);
        end
    end

    numIn  = numel(find_system([mdl '/Electrochem'], 'FindAll','on', ...
                              'LookUnderMasks','all','BlockType','Inport'));
    numOut = numel(find_system([mdl '/Electrochem'], 'FindAll','on', ...
                              'LookUnderMasks','all','BlockType','Outport'));
    fprintf('    Electrochem ports: %d in / %d out\n', numIn, numOut);
    if numIn < 5 || numOut < 3
        close_system(mdl, 0);
        error('MATLAB Function block has %d in / %d out ports (need 5 / 3).', ...
              numIn, numOut);
    end

    % ---- PreLoadFcn defaults (for FUTURE fresh loads) ----
    cb = [ ...
        'if ~exist(''SIM_A'',''var''), SIM_A = 100; end;' newline ...
        'if ~exist(''SIM_N'',''var''), SIM_N =  20; end;' newline ...
        'if ~exist(''SIM_J'',''var''), SIM_J =   1; end;' newline ...
        'if ~exist(''SIM_T'',''var''), SIM_T =  70; end;' newline ...
        'if ~exist(''SIM_P'',''var''), SIM_P =  10; end' ];
    set_param(mdl, 'PreLoadFcn', cb);

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
%SETMATLABFUNCTIONSCRIPT  Populate the MATLAB Function block's script.
%   Tries Stateflow EMChart first (the path that works on R2023b+).
    fullPath = [mdl '/' blockPath];
    ok = false;

    % --- Strategy 1: sfroot EMChart by path ---
    try
        rt = sfroot;
        charts = rt.find('-isa', 'Stateflow.EMChart');
        for i = 1:numel(charts)
            pNorm = strrep(charts(i).Path, '\', '/');
            if strcmpi(pNorm, fullPath) || endsWith(pNorm, ['/' blockPath])
                charts(i).Script = scriptText;
                fprintf('    [setMFScript] EMChart strategy OK (chart %d)\n', i);
                ok = true;
                return;
            end
        end
    catch err1
        fprintf('    [setMFScript] EMChart strategy failed: %s\n', err1.message);
    end

    % --- Strategy 2: MATLABFunctionConfiguration .Script ---
    try
        cfg = get_param(fullPath, 'MATLABFunctionConfiguration');
        cfg.Script = scriptText;
        fprintf('    [setMFScript] MATLABFunctionConfiguration strategy OK\n');
        ok = true;
        return;
    catch err2
        fprintf('    [setMFScript] MFC strategy failed: %s\n', err2.message);
    end

    % --- Strategy 3: set_param on config object ---
    try
        cfg = get_param(fullPath, 'MATLABFunctionConfiguration');
        set_param(cfg, 'Script', scriptText);
        fprintf('    [setMFScript] set_param(config) strategy OK\n');
        ok = true;
        return;
    catch err3
        fprintf('    [setMFScript] set_param strategy failed: %s\n', err3.message);
    end

    % --- Strategy 4: close/reopen then retry EMChart ---
    try
        save_system(mdl);
        close_system(mdl, 0);
        load_system(mdl);
        rt = sfroot;
        charts = rt.find('-isa', 'Stateflow.EMChart');
        for i = 1:numel(charts)
            pNorm = strrep(charts(i).Path, '\', '/');
            if endsWith(pNorm, ['/' blockPath])
                charts(i).Script = scriptText;
                fprintf('    [setMFScript] reopen+EMChart strategy OK (chart %d)\n', i);
                ok = true;
                return;
            end
        end
    catch err4
        fprintf('    [setMFScript] reopen strategy failed: %s\n', err4.message);
    end
end
