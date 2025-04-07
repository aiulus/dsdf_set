function [u_opt, y_opt] = optDDSF(lookup, u_l, traj_ini)
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
        'lambda', 0.01, ...  % Pre-conditioning of the Hankel matrix
        'epsilon', 0.001, ... % For constraint relaxation
        'ridge', 0.01 ...   % Ridge-penalty on norm(alpha, 2)
        );

    %% Extract parameters
    opt_params = lookup.opt_params;

    % Lengths and dimensions
    T_ini = lookup.config.T_ini;
    N = lookup.config.N;    
    L = N + 2 * T_ini;
    m = lookup.dims.m;
    p = lookup.dims.p;
    num_cols = lookup.dims.hankel_cols;
    
    % Matrices
    R = lookup.opt_params.R;
    H = lookup.H;
    T = eye(L, L);
    
    % Initial trajectory 
    u_ini = traj_ini(1:m, :);
    y_ini = traj_ini(m+1:end, :);

    u_e = lookup.sys.S_f.symbolic_solution.u_e;
    y_e = lookup.sys.S_f.symbolic_solution.y_e;

    u_eq = vectorFilter(u_ini(:, end), u_e);
    y_eq = vectorFilter(y_ini(:, end), y_e);

    if iscell(u_eq), u_eq = cell2mat(u_eq); end
    if iscell(y_eq), y_eq = cell2mat(y_eq); end

    %% TODO: S_f can be used to encode a desired target state

    %% Define symbolic variables
    alpha = sdpvar(num_cols, 1);
    control_u = sdpvar(m, L);
    control_y = sdpvar(p, L);
    
    %% Populate the initial and terminal parts of the trajectory
    % Replaces the encodings in constraints
    control_u(:, 1:T_ini) = u_ini;
    control_y(:, 1:T_ini) = y_ini;
    control_u(:, end - T_ini + 1 : end) = repmat(u_eq, 1, T_ini); % T_ini + N + 1 : end
    control_y(:, end - T_ini + 1 : end) = repmat(y_eq, 1, T_ini);

    %% Flatten the variables 
    u_bar = reshape(control_u.', [], 1);
    y_bar = reshape(control_y.', [], 1);
    traj_p_bar = [u_bar; y_bar];

    %% Define the objective function and the constraints  
    delta_u = reshape(control_u(:, 1+T_ini:end-T_ini) - u_l, [], 1);
    objective = delta_u.' * kron(eye(N), R) * delta_u;

    if opt_params.regularize
        switch reg_params.reg_mode
            case 'fro'
                H = low_rank_appr(H);
                H = H / norm(H, 'fro'); % Normalize with Frobenius norm
                lambda = reg_params.lambda;
                H = H + lambda * eye(size(H)); 
            case 'lra' % DOESN'T PRESERVE RANK!!
                H = low_rank_appr(H);
            otherwise
                error('Unknown regularization type: %s', reg_params.reg_mode);
        end
        ridge = reg_params.ridge;
        objective = objective + ridge * (alpha.' * alpha);
    end

    if lookup.opt_params.target_penalty
        target = lookup.sys.params.target;
        target(isnan(target)) = 0;
        target = repmat(target, L, 1);        
        delta_y = y_bar - target;
        delta_y = reshape(delta_y, p, L);

        discount = 0.9;
        gamma = discount .^ (L:-1:1);
        gamma = repmat(gamma, p, 1);
        delta_y = delta_y .* gamma;

        tpt = trace(delta_y * T * delta_y.');
        objective = objective + tpt;
    end


    if ~opt_params.regularize
        epsilon = 0; 
    else 
        epsilon = reg_params.epsilon; 
    end
    
    constraints = traj_p_bar == H * alpha;

    U_set = lookup.sys.sets.constraints.U;
    Y_set = lookup.sys.sets.constraints.Y;

    for k = T_ini+1 : L-T_ini
        switch lookup.opt_params.constr_type
            case 'u'
                constraints = [constraints, zonoContains(U_set, control_u(:, k))];
            case 'y'
                constraints = [constraints, zonoContains(Y_set, control_y(:, k))];
            case 'f'
                constraints = [constraints, ...
                    zonoContains(U_set, control_u(:, k)), ...
                    zonoContains(Y_set, control_y(:, k))];
        end
    end

    % === Terminal constraint)
    y_terminal = control_y(:, L - T_ini + 1);
    constraints = [constraints, zonoContains(lookup.sys.sets.S_f_zono, y_terminal)];

    % === Solve QP
    switch opt_params.solver_type
        case 'q'
            options = sdpsettings('solver', 'quadprog', 'verbose', 0);
        case 'f'
            options = sdpsettings('solver', 'fmincon', 'verbose', 0);
        otherwise
            options = sdpsettings('solver', 'OSQP', ...
                'osqp.max_iter', 20000, ...
                'osqp.eps_abs', 1e-3, ...
                'osqp.eps_rel', 1e-3, ...
                'osqp.polish_refine_iter', 2, ...
                'osqp.scaling', 10, ...
                'warmstart', 0, ...
                'verbose', lookup.IO_params.verbose);
    end


    diagnostics = optimize(constraints, objective, options);
        
    if diagnostics.problem == 0 % Feasible solution found
        % Extract optimal values
        u_opt = value(control_u);
        y_opt = value(control_y);
    else
        % fprintf('\n---------------------------- OPTIMIZER FAILED ----------------------------\n');
        disp(diagnostics.problem); % Solver exit code
        disp(diagnostics.info);    % Detailed solver feedback        
    end

    if lookup.IO_params.verbose
        disp('---- Debug: Objective Function Evaluation ----');
        disp('Objective value:');
        disp(value(objective)); % Ensure it computes as expected
        disp('---- Debug: Feasibility Check ----');
        disp('Max constraint violation:');
        disp(max(check(constraints))); % Show largest constraint violation
    end

end


function C = zonoContains(Z, x)
    % Zonotope Z: Z = [c | G]
    G = Z.Z(:, 2:end);  % Generators
    c = Z.Z(:, 1);      % Center
    n_gen = size(G, 2);
    alpha = sdpvar(n_gen, 1);
    C = [G * alpha + c == x, norm(alpha, Inf) <= 1];
end


