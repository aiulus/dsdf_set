function [lookup, time, logs] = runDDSF(systype, T_sim, N, T_ini, scale_constraints, R, toggle_plot)
    toggle_save = true;

    if nargin < 7
        toggle_plot = 0;
    end

    % Input signal generation options:
    % {'prbs', 'sinusoid', 'sinusoidal_sweep', 'uniform', 
    % 'custom_uniform', 'controlled_random', 'white_noise'}

    %% Step 1: Configuration
    data_options = struct( ...
        'datagen_mode', 'controlled_random', ... 
        'scale', 3, ... % Constraint relaxation for the learning agent
        'safe', false ... % Set 1/true to ensure the safety of u_d
        );
    
    run_options = struct( ...
        'system_type', systype, ...
        'T_sim', T_sim, ...
        'T_d', 0 ... % Input delay [s] / dt [s]
        );
    
    IO_params = struct( ...
        'debug', true, ...
        'save', true, ...
        'log_interval', 1, ... % Log every n-the step of the simulation
        'verbose', false ...
        );
    
    % TODO: Add information on various the configuration options
    opt_params = struct( ...                       
                        'regularize', true, ...
                        'constr_type', 'f', ...
                        'solver_type', 'o', ...
                        'target_penalty', false, ...    
                        'init', true, ... % Encode initial condition
                        'R', 1 ...
                       );
    
    if ismember(systype, {'test_nonlinear', 'van_der_pol', 'nonlinear_pendulum'})
        sys = nonlinearSysInit(systype);
    else
        sys = systemsDDSF(systype);
    end

    dims = sys.dims;

    % Reasign safety boundaries if specidied in the input
    if scale_constraints ~= -1
        sys.constraints.U = updateBounds(sys.constraints.U, scale_constraints);
        sys.constraints.Y = updateBounds(sys.constraints.Y, scale_constraints);
    end

    % Match T_ini, N to default values if not specified in the input
    if T_ini == -1 || N == -1
        T_ini = sys.config.T_ini;
        N = sys.config.N;
    end

    % Reassign weight matrix if specified in the input
    if R ~= -1
        opt_params.R = R;
    end
    % Upscale to correct dimensions
    opt_params.R = opt_params.R * eye(dims.m);
        
    % Create struct object 'lookup' for central and extensive parameter passing.
    lookup = struct( ...
                    'sys', sys, ...
                    'systype', systype, ...
                    'opt_params', opt_params, ...
                    'config', sys.config, ...
                    'dims', dims, ...
                    'IO_params', IO_params, ...
                    'T_sim', run_options.T_sim, ...
                    'data_options', data_options, ...
                    'T_d', run_options.T_d ...
                    );

    lookup.config.N = N;
    lookup.config.T_ini = T_ini;    

    lookup.sys.config.N = N;
    lookup.sys.config.T_ini = T_ini;    

    %% Step 2: Generate data & Hankel matrices
    [u_d, y_d, ~, ~, ~] = gendataDDSF(lookup); 
    
    [H_u, H_y] = hankelDDSF(u_d, y_d, lookup);
    
    %% TODO: H_v_2 = construct_hankel(v_d, T_ini + 1)
    H_u_2 = construct_hankel(u_d, 2);
    H_y_2 = construct_hankel(y_d, 2);
     
    lookup.H = [H_u; H_y];
    lookup.H_u = H_u; lookup.H_y = H_y;
    lookup.dims.hankel_cols = size(H_u, 2);
        
    %% Initialize objects to log simulation history
    T_ini = lookup.config.T_ini;
    logs = struct( ...
            'u', [u_d(:, 1:T_ini).'; ...
            zeros(dims.m, T_sim).'].', ... 
            'y', [y_d(:, 1:T_ini).'; ...
            zeros(dims.p, T_sim).'].', ... 
            'yl', zeros(dims.p, T_sim), ...
            'ul', zeros(dims.m, lookup.config.N, T_sim), ...
            'ul_t', zeros(dims.m, T_sim), ...
            'loss', zeros(2, T_ini + T_sim) ...
        );
    
    % Obtain an L-step random control policy
    T_d = run_options.T_d;
    %% Step 3: Receding Horizon Loop
    for t=(T_ini + 1):(T_ini + T_sim)
        fprintf("----------------- DEBUG Information -----------------\n");
        fprintf("CURRENT SIMULATION STEP: t = %d\n", t - T_ini);
    
        u_l = learning_policy(lookup);
        logs.ul(:, :, t - T_ini) = u_l;
        ul_t = u_l(:, 1);
        logs.ul_t(:, t - T_ini) = ul_t;
    
        u_ini = logs.u(:, (t - T_ini):(t-1));
        y_ini = logs.y(:, (t - T_ini):(t-1));
        traj_ini = [u_ini; y_ini]; 
    
        [u_opt, y_opt] = optDDSF(lookup, u_l, traj_ini);
        
        % loss_t = get_loss(lookup, ul_t, u_opt, y_opt);
    
        u_next = u_opt(:, 1 + T_ini);
        y_next = y_opt(:, 1 + T_ini);
                
        logs.ul(:, :, t) = u_l;
        logs.u(:, t) = u_next;
        logs.y(:, t) = y_next;
        % logs.loss(:, t) = loss_t;


        %% TODO: if (t - T_ini + 1) >= 1 && T_ini ~= 1, try ... catch ... end; end
        %%      yl_next = dataBasedU2Y(ul_t, logs.u(:, t- T_ini + 1 : t), logs.y(:, t- T_ini + 1 : t), H_u_2, H_y_2);
        try            
            yl_next = dataBasedU2Y(ul_t, logs.u(:, t),logs.y(:, t), H_u_2, H_y_2);
            logs.yl(:, t-T_ini) = yl_next;
        catch ME
            fprintf('Failed to execute dataBasedU2Y at time step %d: %s', t-T_ini, ME.message);
        end
    end
    
    % Store the final simulation results
    logs.u = logs.u(:, 1+T_ini:end);
    logs.y = logs.y(:, 1+T_ini:end);
    lookup.logs = logs;
    
    %% Plot the results
    time = 1:T_sim;
    if toggle_plot
        plotDDSF(time, logs, lookup)
    end 
    if toggle_save
        prefix_u = sprintf(strcat('U-ddsf-', systype, '-single_run'));
        prefix_y = sprintf(strcat('Y-ddsf-', systype, '-single_run'));
        output_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'outputs', 'data');
        if ~exist(output_dir, 'dir'), mkdir(output_dir); end
        fullpath_u = fullfile(output_dir, prefix_u);
        fullpath_y = fullfile(output_dir, prefix_y);
        csvFlexSave(fullpath_u, logs.u, logs.ul_t);
        csvFlexSave(fullpath_y, logs.y, logs.yl);
    end
end

function loss = get_loss(lookup, u_l, u_opt, y_opt)
    R = lookup.opt_params.R;
    y_target = lookup.sys.params.target;

    loss1 = det(R) * norm(u_l - u_opt);    
    loss2 = loss1 + norm(y_target - y_opt);
    
    loss = [loss1; loss2];
end