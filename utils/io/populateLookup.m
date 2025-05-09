function lookup = populateLookup(systype, T_sim, N, T_ini, scale_constraints, R, toggle_plot, x0_mode)
    if nargin < 8, x0_mode = 'single'; end

    % === Configuration ===
    data_options = struct( ...
        'datagen_mode', 'controlled_random', ...
        'safe', false ...
    );

    opt_params = struct( ...
        'regularize', true, ...
        'constr_type', 'f', ...
        'solver_type', 'o', ...
        'target_penalty', false, ...
        'init', true, ...
        'R', 1 ...
    );

    uncertainty_options = struct('x0_mode', x0_mode);

    if ismember(systype, {'test_nonlinear', 'van_der_pol', 'nonlinear_pendulum'})
        sys = nonlinearSysInit(systype);
    else
        sys = systemsDDSF(systype);
    end

    if scale_constraints ~= -1
        sys.constraints.U = updateBounds(sys.constraints.U, scale_constraints);
        sys.constraints.Y = updateBounds(sys.constraints.Y, scale_constraints);
        sys = system2set(sys);
    end

    if T_ini == -1 || N == -1
        T_ini = sys.config.T_ini;
        N = sys.config.N;
    end

    dims = sys.dims;
    if R ~= -1
        opt_params.R = R;
    end
    opt_params.R = opt_params.R * eye(dims.m);

    IO_params = struct( ...
        'debug', true, ...
        'save', false, ...
        'log_interval', 1, ...
        'verbose', false ...
    );

    lookup = struct( ...
        'sys', sys, ...
        'systype', systype, ...
        'opt_params', opt_params, ...
        'config', sys.config, ...
        'dims', dims, ...
        'IO_params', IO_params, ...
        'T_sim', T_sim, ...
        'data_options', data_options, ...
        'T_d', 0, ...
        'uncertainty_options', uncertainty_options ...
    );

    lookup.config.N = N;
    lookup.config.T_ini = T_ini;
    lookup.sys.config.N = N;
    lookup.sys.config.T_ini = T_ini;
end
