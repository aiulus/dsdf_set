function [u, y] = deepc3Opt(lookup, H, u_ini, y_ini)
    % reg_params.reg_mode   -   OPTIONS
    %                       - 'fro': Uses Frobenius norm + Tikhonov
    %                                regularization on the Hankel matrix
    %                       - 'svd': Uses SVD-low rank approximation on the
    %                                Hankel matrix
    % reg_params.lambda     -  Tikhonov regularization parameter -
    %                          Only to be used in mode 'fro'
    % reg_params.epsilon    -  For constraint relaxation
    % reg_params.ridge      -  Penalizes higher norm-values of alpha
    
    reg_params = struct( ...
        'reg_mode', 'fro', ... % Regularization mode
        'lambda', 0.01, ...
        'epsilon', 0.01, ...
        'ridge', 0.01 ...
        );    
    
    verbose = false; % Toggle debug mode
    optimizer_type = 'o'; % Toggle optimization type
    constr_type = 'f'; % Toggle constraint type
    
    %% Extract DeePC parameters
    % TODO: Line 8 gets interpreted as a function call by the b&b solver
    Q = lookup.deepc.Q;
    R = lookup.deepc.R;
    
    % Extract dimensions
    m = lookup.dims.m;
    p = lookup.dims.p;
    
    % Extract configuration parameters
    T = lookup.config.T;
    T_ini = lookup.config.T_ini;
    N = lookup.config.N;
    
    % Extract (past/future slices of) Hankel matrices
    Up = H.Up; Yp = H.Yp;
    Uf = H.Uf; Yf = H.Yf;
    
    % Extract system constraints
    U = lookup.sys.constraints.U;
    Y = lookup.sys.constraints.Y;
    target = lookup.sys.params.target;
    
    %% Define symbolic variables
    g = sdpvar(T - T_ini - N + 1, 1);
    u = sdpvar(N * m, 1);
    y = sdpvar(N * p, 1);    
    
    %% Construct the reference trajectory
    target = reshape(repmat(target, N, 1).', [], 1);
    % nancols = isnan(target);
    numcols = ~isnan(target);
    
    %% Define the optimization objective
    delta = y(numcols) - target(numcols);
    Qext = kron(eye(N), Q);
    Qext_cut = Qext(numcols, numcols);
    
    objective = delta' * Qext_cut * delta + u' * kron(eye(N), R) * u;

    if lookup.opt_params.regularize
        % H = regHankelDDSF(lookup.H_u, lookup.H_y);
        % num_cols = size(H, 2);
        switch reg_params.reg_mode
            case 'fro'
                H = H / norm(H, 'fro'); % Normalize with Frobenius norm
                lambda = reg_params.lambda;
                H = H + lambda * eye(size(H)); % Tikhonov regularization
            case 'svd'
                H = low_rank_appr(H);
            otherwise
                error('Unknown regularization type: %s', reg_params.reg_mode);
        end

        ridge = reg_params.ridge;
        objective = objective + ridge * (alpha.' * alpha);
    end
    
    u_min = U(:, 1); u_max = U(:, 2);
    y_min = Y(:, 1); y_max = Y(:, 2);    
    u_min(isinf(u_min)) = 1e+8; u_max(isinf(u_max)) = 1e+8;
    y_min(isinf(y_min)) = 1e+8; y_max(isinf(y_max)) = 1e+8;

    if ~lookup.opt_params.regularize
        epsilon = 0; 
    else 
        epsilon = reg_params.epsilon; 
    end
    
    epsilon = epsilon / 2;
    delta_u = epsilon .* (u_max - u_min);
    delta_y = epsilon .* (y_max - y_min);
    
    u_lb = reshape(repmat(u_min - delta_u, N, 1).', [], 1);
    u_ub = reshape(repmat(u_max + delta_u, N, 1).', [], 1);
    y_lb = reshape(repmat(y_min - delta_y, N, 1).', [], 1);
    y_ub = reshape(repmat(y_max + delta_y, N, 1).', [], 1);
    
    % Define the constraints
    if lower(constr_type) == 'f'           % Full
        constraints = [Up * g == u_ini, ...
            Yp * g == y_ini, ...
            Uf * g == u, ...
            Yf * g == y, ...
            u_lb <= u, u <= u_ub, ...
            y_lb <= y, y <= y_ub ...
            ];
    elseif lower(constr_type) == 's'       % Simple, just models the system
        constraints = [Up * g == u_ini, ...% dynamics - same as setting
            Yp * g == y_ini, ...% all input constraints to
            Uf * g == u, ...    % +/-inf
            Yf * g == y
            ];
    elseif lower(constr_type) == 'e'       % Empty - meant to be used just
        constraints = [];                  % as sanity check
    end
    
    if lower(optimizer_type) == 'q'
        options = sdpsettings('solver', 'quadprog', 'verbose', verbose);
    elseif lower(optimizer_type) == 'f'
        options = sdpsettings('solver', 'fmincon', ...
            'sdpa.maxIteration', 10, ...
            'verbose', 0);
    elseif lower(optimizer_type) == 'o'
        %options = sdpsettings('solver', 'OSQP', ... % Use OSQP solver
        %    'verbose', verbose, ...
        %    'osqp.max_iter', 30000, ...   % Set maximum iterations
        %    'osqp.eps_abs', 1e-5, ...     % Absolute tolerance
        %    'warmstart', 0);             % Disable warm start
        options = sdpsettings('solver', 'OSQP', ...
                  'verbose', verbose, ...             % Detailed solver output
                  'osqp.max_iter', 30000, ...         % Set maximum iterations
                  'osqp.eps_abs', 1e-5, ...           % Absolute tolerance
                  'osqp.eps_rel', 1e-5, ...           % Relative tolerance
                  'osqp.adaptive_rho', true, ...      % Enable adaptive rho
                  'osqp.adaptive_rho_interval', 50, ...
                  'osqp.adaptive_rho_tolerance', 1, ...
                  'osqp.polish_refine_iter', 2, ...
                  'osqp.scaling', 10, ...             % Number of scaling iterations
                  'warmstart', 0);                    % Disable warm start
    elseif lower(optimizer_type) == 'b'
        options = sdpsettings('solver', 'bmibnb', 'verbose', 1);
    else
        error('Error assigning solver options!');
    end
    
    diagnostics = optimize(constraints, objective, options);
    
    if diagnostics.problem == 0
        % g = value(g);
        u = value(u);
        y = value(y);
        if verbose
            disp("REFERENCE TRAJECTORY: "); disp(target');
            disp("DELTA: "); disp(value(delta)');
            disp("PEANALTY TERM: "); disp(value(objective));
        end
    else
        disp(diagnostics.info);
        error('The problem did not solve successfully.');
    end
end
    
    %   CRUISE CONTROL - DeePC
    %
    %   best_params =
    %
    %  struct with fields:
    %
    %       Q: 10000
    %        R: 0.0100
    %    T_ini: 10
    %       N: 30
    %
    %   in 1348/60