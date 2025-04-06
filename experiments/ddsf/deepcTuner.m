%% TODO: Remove
function [u, y, descriptions] = deepcTuner(mode, systype, T_sim, toggle_save)
    % RETURNS -
    %       u -         R^(n_runs x m)
    %       y -         R^(n_runs x p)
    % descriptions -    n_runs x strings naming the respective parameter
    %                   configuration

    max_tries = 5;
    
    if nargin < 4
        toggle_save = false;
    end
    
    Qs = logspace(-2, 4, 7);
    Rs = logspace(-2, 2, 5);
    
    vals = struct( ...
        'NvsTini', [ ...
        1 * ones(6, 1), (5:5:30)'; ...
        2 * ones(6, 1), (5:5:30)'; ...
        5 * ones(5, 1), (10:5:30)'; ...
        10 * ones(4, 1), (15:5:30)' ...
        ], ...
        'QvsR', table2array(combinations(Qs, Rs)), ...
        'mixed', struct( ...
        'nt', [ ...
        2 * ones(6, 1), (5:5:30)'; ...
        5 * ones(5, 1), (10:5:30)' ...
        ], ...
        'qr', table2array(combinations(logspace(-2, 1, 4), logspace(-2, 1, 4))) ...
        ) ...
        );
    
    
    switch mode
        case {'NvsTini', 'nt'}
            values = vals.NvsTini;
            nruns = max(size(values));
        case {'QvsR', 'qr'}
            values = vals.QvsR;
            nruns = max(size(values));
        case 'mixed'
            values = vals.mixed;
            nruns = max(size(values.nt)) * max(size(values.qr));
        otherwise
            error("The provided tuner mode isn't supported.");
    end
    
    u = cell(1, nruns); y = cell(1, nruns); descriptions = cell(1, nruns);
    
    output_dir = fullfile('..', 'outputs', 'plots', 'deepc_tuner');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    switch mode
        case {'NvsTini', 'nt'}
            Q = -1; R = -1;
            tic;
            for i=1:nruns
                fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);
                if i > 1
                    elapsed = toc;
                    timeEstimator(elapsed, i, nruns);
                end
                T_ini_i = values(i, 1);
                N_i = values(i, 2);
                d_i = sprintf('deepcTuner-%s-N%d-Tini%d-%s-T%d', mode, N_i, T_ini_i, systype, T_sim);
                for k=1:max_tries
                    try
                        logs = runParamDPC(systype, Q, R, T_ini_i, N_i, T_sim);
                        u_i = logs.u_sim; % m x T_sim
                        y_i = logs.y_sim; % m x T_sim
                        u{i} = u_i;
                        y{i} = y_i;
                        descriptions{i} = d_i;
                        break;
                    catch ME
                        fprintf(['Attempt to run runParamDPC.m (conf.: %s) failed at: %s\n ' ...
                            'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                    end
                end
            end
        case {'QvsR', 'qr'}
            T_ini = -1; N = -1;
            tic;
            for i=1:nruns
                fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);
                if i > 1
                    elapsed = toc;
                    timeEstimator(elapsed, i, nruns);
                end                
                Q_i = values(i, 1);
                R_i = values(i, 2);
                d_i = sprintf('deepcTuner-%s-Q%d-R%d-%s-T%d', mode, Q_i, R_i, systype, T_sim);
                for k=1:max_tries
                    try
                        logs = runParamDPC(systype, Q_i, R_i, T_ini, N, T_sim);
                        u_i = logs.u_sim; % m x T_sim
                        y_i = logs.y_sim; % m x T_sim
                        u{i} = u_i;
                        y{i} = y_i;
                        descriptions{i} = d_i;
                        break;
                    catch ME
                        fprintf(['Attempt to run DDSF (conf.: %s) failed at: %s\n ' ...
                            'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                    end
                end
            end
        case 'mixed'
            qr = values.qr;
            nt = values.nt;
            for i=1:max(size(qr))
                Q_i = qr(i, 1); R_i = qr(i, 2);
                tic;
                for j=1:max(size(nt))                    
                    T_ini_j = nt(j, 1); N_j = nt(j, 2);
                    d_i = sprintf('deepcTuner-%s-N%d-Tini%d-Q%d-R%d-%s-T%d', mode, N_j, T_ini_j, Q_i, R_i, systype, T_sim);
                    t = (i-1)*max(size(nt))+j;
                    fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', t, nruns);
                    if t > 1
                        elapsed = toc;
                        timeEstimator(elapsed, t, nruns);
                    end
                    for k=1:max_tries
                        try
                            logs = runParamDPC(systype, Q_i, R_i, T_ini_j, N_j, T_sim);
                            u_i = logs.u_sim; % m x T_sim
                            y_i = logs.y_sim; % m x T_sim
                            u{i} = u_i;
                            y{i} = y_i;
                            descriptions{i} = d_i;
                            break;
                        catch ME
                            fprintf(['Attempt to run runParamDPC.m (conf.: %s) failed at: %s\n ' ...
                                'Message: = %s\n. Trying again...'], d_i, ME.stack(1).name, ME.message);
                        end
                    end
                end
            end
    end

    if toggle_save
        prefix = sprintf('deepcTuner-%s-%s-T%d', mode, systype, T_sim);
        % save2csv(u, y, descriptions, prefix);
        csvFlexSave(prefix, u,  y, descriptions);
    end

    u = cell2mat(u(:));
    y = cell2mat(y(:));
end

