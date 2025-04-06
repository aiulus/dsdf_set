function [time, u_sim, y_sim] = runMPC(systype, T_sim)
    %% Define the system
    sys = sysInit(systype);
    verbose = true; 
    save = true;
    rng(0, 'twister'); % Set seed and generator
    
    %% TODO: Ensure shape match across _systems & Opt !!
    
    %% Extract relevant parameters
    dims = sys.dims;
    config =  sys.config;
    config.T_sim = T_sim;
    s = 1;
    config.s = s;
    opt_params = sys.opt_params;
    
    %% Initialize containers for logging
    u_sim = zeros(dims.m, config.T_sim + 1);
    y_sim = zeros(dims.p, config.T_sim + 1);
    x_sim = zeros(dims.n, config.T_sim + 1);
    
    
    %% Data structure that holds all relevant parameters
    lookup = struct( ...
            'dims', dims, ...
            'config', config, ...
            'opt_params', opt_params, ...
            'sys', sys, ...
            'T_d', 0 ... % Input delay [s] / dt [s]
            );
    
    x_ini = sys.params.x_ini;
    x_sim(:, 1) = x_ini;
    x = x_ini;
    
    %% Receding Horizon Loop
    for t=1:config.T_sim
        % Solve the quadratic optimization problem
        [u_opt, ~, ~] =  mpc3Opt(lookup, x);
    
        u = u_opt(1:dims.m*config.s);
        %y = y_opt(1:dims.p*config.s);
        %x = x_opt(1:dims.n*config.s);
        x = sys.A*x + sys.B*u;
        y = sys.C*x + sys.D*u;
    
        % Log the resulting trajectory
        u_sim(:, t) = u;
        y_sim(:, t) = y;
        x_sim(:, t + 1) = x;
    
        fprintf("------------------- " + ...
                 "MPC - Simulation step %d / %d -------------------\n", t, config.T_sim);
        if verbose
            fprintf("Distance to target: "); disp(norm(y - sys.params.target));
        end
    end
    
    %% Plot the results
    Tsim = config.T_sim; 
    time = 0:Tsim;
    %deepc3plot(time, y_sim, u_sim, sys)
    
    %% Save the results
    if save
        save2csv(time, u_sim, y_sim);
    end

end