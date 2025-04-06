function [lookup, time, logs] = singleRunDDSF(systype, T_sim, toggle_plot)
    % Input signal generation options:
    % {'prbs', 'sinusoid', 'sinusoidal_sweep', 'uniform', 
    % 'custom_uniform', 'controlled_random', 'white_noise'}

    %% Step 1: Configuration
    data_options = struct( ...
        'datagen_mode', 'controlled_random', ...
        'scale', 1.5, ... % Constraint relaxation for the learning agent
        'safe', false ... % Set 1/true to sample from safe input set
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
    opt_params.R = opt_params.R * eye(dims.m);    
    
    % Create struct object 'lookup' for central and extensive parameter passing.
    lookup = struct( ...
                    'systype', systype, ...
                    'sys', sys, ...
                    'opt_params', opt_params, ...
                    'config', sys.config, ...
                    'dims', dims, ...
                    'IO_params', IO_params, ...
                    'T_sim', run_options.T_sim, ...
                    'data_options', data_options, ...
                    'T_d', run_options.T_d ...
                    );  

    %% Step 2: Generate data & Hankel matrices
    [u_d, y_d, ~, ~, ~] = gendataDDSF(lookup); 
    
    [H_u, H_y] = hankelDDSF(u_d, y_d, lookup);

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
        
        loss_t = get_loss(lookup, ul_t, u_opt, y_opt);
    
        u_next = u_opt(:, 1 + T_ini);
        y_next = y_opt(:, 1 + T_ini);
                
        logs.ul(:, :, t) = u_l;
        logs.u(:, t) = u_next;
        logs.y(:, t) = y_next;
        logs.loss(:, t) = loss_t;

        try            
            yl_next = dataBasedU2Y(ul_t, logs.u(:, t-1),logs.y(:, t-1), H_u_2, H_y_2);
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

end

function loss = get_loss(lookup, u_l, u_opt, y_opt)
    R = lookup.opt_params.R;
    y_target = lookup.sys.params.target;

    loss1 = det(R) * norm(u_l - u_opt);    
    loss2 = loss1 + norm(y_target - y_opt);
    
    loss = [loss1; loss2];
end