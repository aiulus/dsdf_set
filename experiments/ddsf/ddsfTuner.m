%% TODO: Remove
function [u, ul, y, yl, descriptions] = ddsfTuner(mode, systype, T_sim, toggle_save)
    % RETURNS -
    %       u -         R^(n_runs x m)
    %       ul -        R^(n_runs x m)
    % descriptions -    n_runs x strings naming the respective parameter
    %                   configuration
    
    max_tries = 2;
    toggle_plot = false;
    
    if nargin < 4
        toggle_save = false;
    end
    

    vals = struct( ...
    'r', 10.^(-8:1:8), ...
    'NvsTini', [ ...
    1 * ones(6, 1), (5:5:30)'; ...
    2 * ones(6, 1), (5:5:30)'; ...
    5 * ones(5, 1), (10:5:30)'; ...
    10 * ones(4, 1), (15:5:30)' ...
    ], ...
    'constraints', [1e+8, 0.1, 0.5, 1.5, 2, 10, 100], ...
    'mixed', struct( ...
    'nt', [vertcat(repelem([2; 3; 5], 4), 10), vertcat(repmat(5:5:20, 1, 2)', (10:5:25).', 30)], ...
    'constr', [0.1, 0.5, 1, 1.5, 2, 10, 1e+8], ...
    'R', [1e-4, 0.1, 1, 10, 100, 1e+6] ...
    ) ...
    );
    
    switch lower(mode)
        case 'r'
            values = vals.r;
            nruns = max(size(values));
        case {'nvstini', 'nt'}
            values = vals.NvsTini;
            nruns = max(size(values));
        case {'constraints', 'constr'}
            values = vals.constraints;
            nruns = max(size(values));
        case 'mixed'
            values = vals.mixed;
            nruns = max(size(values.nt)) * max(size(values.constr)) * max(size(values.R));
        otherwise
            error("The provided tuner mode isn't supported.");
    end
    
    u = cell(1, nruns); ul = cell(1, nruns); descriptions = cell(1, nruns);
    y = cell(1, nruns); yl = cell(1, nruns);
    
    output_dir = fullfile('..', 'outputs', 'plots', 'ddsf_tuner');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    switch mode
        case 'r'
            constraint_scaler = 2; T_ini = -1; N = -1;
            tic;
            for i=1:nruns
                R_i = values(i);
                d_i = sprintf('ddsfTuner-%s-N%d-Tini-%d-R%d-%s-T%d', mode, N, T_ini, R_i, systype, T_sim);
                fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);
                if i > 1
                    elapsed = toc;
                    timeEstimator(elapsed, i, nruns);
                end
                for k=1:max_tries
                    fprintf(['------------------- Attempt %d / %d for configuration %d / %d ' ...
                        '-------------------\n'], k, max_tries, i, nruns);
                    try
                        [~, ~, logs] = runDDSF(systype, T_sim, N, T_ini, constraint_scaler, R_i, toggle_plot);
    
                        u_i = logs.u; % m x T_sim
                        ul_i = logs.ul_t; % m x T_sim
    
                        u{i} = u_i;
                        ul{i} = ul_i;
    
                        y_i = logs.y;
                        yl_i = logs.yl;
    
                        y{i} = y_i;
                        yl{i} = yl_i;
    
                        descriptions{i} = d_i;
                        break;
                    catch ME
                        fprintf(['Attempt to run runDDSF.m (conf.: %s) failed at: %s\n ' ...
                            'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                    end
                end
            end
        case {'NvsTini', 'nt'}
            constraint_scaler = 1; R = -1;
            tic;
            for i=1:nruns
                T_ini_i = values(i, 1);
                N_i = values(i, 2);
                d_i = sprintf('ddsfTuner-%s-N%d-Tini-%d-%s-T%d', mode, N_i, T_ini_i, systype, T_sim);
                fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);
                if i > 1
                    elapsed = toc;
                    %timeEstimator(elapsed, i, nruns);
                end
                for k=1:max_tries
                    try
                        [~, ~, logs] = runDDSF(systype, T_sim, N_i, T_ini_i, constraint_scaler, R, toggle_plot);
    
                        u_i = logs.u; % m x T_sim
                        ul_i = logs.ul_t; % m x T_sim
    
                        u{i} = u_i;
                        ul{i} = ul_i;
    
                        y_i = logs.y;
                        yl_i = logs.yl;
    
                        y{i} = y_i;
                        yl{i} = yl_i;
    
                        descriptions{i} = d_i;
                        break;
                    catch ME
                        fprintf(['Attempt to run runDDSF.m (conf.: %s) failed at: %s\n ' ...
                            'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                    end
                end
            end
        case {'constraints', 'constr'}
            T_ini = -1; N = -1;
            tic;
            for i=1:nruns
                constr_i = values(i); R = -1;
                d_i = sprintf('ddsfTuner-%s-scalingfactor%d-%s-T%d', mode, values(i), systype, T_sim);
                fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);
                if i > 1
                    elapsed = toc;
                    timeEstimator(elapsed, i, nruns);
                end
                for k=1:max_tries
                    try
                        [~, ~, logs] = runDDSF(systype, T_sim, N, T_ini, constr_i, R, toggle_plot);
    
                        u_i = logs.u; % m x T_sim
                        ul_i = logs.ul_t; % m x T_sim
    
                        u{i} = u_i;
                        ul{i} = ul_i;
    
                        y_i = logs.y;
                        yl_i = logs.yl;
    
                        y{i} = y_i;
                        yl{i} = yl_i;
    
                        descriptions{i} = d_i;
                    catch ME
                        fprintf(['Attempt to run runDDSF.m (conf.: %s) failed at: %s\n ' ...
                            'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                    end
                end
            end
        case 'mixed'
            nt = values.nt;
            constr = values.constr;
            R = values.R;
            tic;
            for i=1:max(size(nt))
                T_ini_i = nt(i, 1);
                N_i = nt(i, 2);
                for j=1:max(size(constr))
                    constr_j = constr(j);
                    for l=1:max(size(R))
                        R_l = R(l);
                        d_i = sprintf('ddsfTuner-%s-Tini%d-N%d-constr_scale%d-R%d-systype-%s-T%d', mode, T_ini_i, N_i,constr_j, R_l, systype, T_sim);
                        t = (i-1)*max(size(constr))+(j-1)*max(size(R))+l;
                        fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', t, nruns);
                        if t > 1
                            elapsed = toc;
                            %timeEstimator(elapsed, t, nruns);
                        end
                        for k=1:max_tries
                            try
                                [~, ~, logs] = runDDSF(systype, T_sim, N_i, T_ini_i, constr_j, R_l, toggle_plot);
    
                                u_i = logs.u; % m x T_sim
                                ul_i = logs.ul_t; % m x T_sim
    
                                u{i} = u_i;
                                ul{i} = ul_i;
    
                                y_i = logs.y;
                                yl_i = logs.yl;
    
                                y{i} = y_i;
                                yl{i} = yl_i;
    
                                descriptions{i} = d_i;
                                break;
                            catch ME
                                fprintf(['Attempt to run runDDSF.m (conf.: %s) failed at: %s\n ' ...
                                    'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                            end
                        end
                    end
                end
            end
    end
    
    if toggle_save
        prefix_u = sprintf('U-ddsfTuner-systype-%s-mode-%s-T%d', systype, mode, T_sim);
        prefix_y = sprintf('Y-ddsfTuner-systype-%s-mode-%s-T%d', systype, mode, T_sim);
        csvFlexSave(prefix_u, u, ul, descriptions);
        csvFlexSave(prefix_y, y, yl, descriptions);
    end
    
    u = cell2mat(u(:));
    ul = cell2mat(ul(:));
end