function [u, x, y] = mpc3Opt(lookup, x_ini)
    verbose = true; % Toggle debug mode
    optimizer_type = 'o'; % Toggle optimization type 
    constr_type = 'f'; % Toggle constraint type

    %% Extract DeePC parameters
    % TODO: Line 8 gets interpreted as a function call by the b&b solver
    Q = lookup.opt_params.Q;
    R = lookup.opt_params.R;

    % Extract dimensions
    m = lookup.dims.m;
    p = lookup.dims.p;
    n = lookup.dims.n;
    
    N = lookup.config.N;

    % Extract system constraints
    U = lookup.sys.constraints.U;
    Y = lookup.sys.constraints.Y;
    target = lookup.sys.params.target;

    %% Define symbolic variables
    u = sdpvar(N, m);  u = reshape(u, [], 1);
    y = sdpvar(N, p); y = reshape(y, [], 1);
    x = sdpvar(N + 1, n); x = reshape(x, [], 1);
    % x(1, :) = x_ini; x = reshape(x, [], 1);

    %% Construct the reference trajectory
    target = reshape(repmat(target, N, 1).', [], 1);
    % nancols = isnan(target);
    numcols = ~isnan(target);    

    %% Define the optimization objective
    delta = y(numcols) - target(numcols);
    Qext = kron(eye(N), Q);
    Qext_cut = Qext(numcols, numcols);

    A = lookup.sys.A;
    B = lookup.sys.B;
    C = lookup.sys.C;
    D = lookup.sys.D;
    
    objective = delta' * Qext_cut * delta + u' * kron(eye(N), R) * u;
    
    u_lb = reshape(repmat(U(:, 1), N, 1).', [], 1);
    u_ub = reshape(repmat(U(:, 2), N, 1).', [], 1);
    y_lb = reshape(repmat(Y(:, 1), N, 1).', [], 1);
    y_ub = reshape(repmat(Y(:, 2), N, 1).', [], 1);

    constraints = x(1:n) == x_ini;

    for i=1:N
        %fprintf("size of B: (%d, %d)\n", size(B, 1), size(B, 2));
        %fprintf("Size of u: "); disp(size(u((i-1)*m + 1:i*m)));
        constraints = [constraints, ...
                       x(i*n + 1 : (i + 1)*n) == A*x((i-1)*n + 1:i*n) + B*u((i-1)*m + 1:i*m), ...
                       y((i-1)*p + 1:i*p) == C*x((i-1)*n + 1:i*n) + D*u((i-1)*m + 1:i*m)];
    end

    % Define the constraints
    if lower(constr_type) == 'f'           % Full 
        constraints = [constraints, ...
                       u_lb <= u & u <= u_ub, ...
                       y_lb <= y & y <= y_ub ...
                      ];
    elseif lower(constr_type) == 'u'       % Only input constraints
        constraints = [constraints, ...
                       u_lb <= u & u <= u_ub];
    elseif lower(constr_type) == 'y'       % Only output constraints
        constraints = [constraints, ...
                       y_lb <= y & y <= y_ub ...
                      ];
    elseif lower(constr_type) == 's'       % Simple, just models the system 
        % Nothing to do
    end

    if lower(optimizer_type) == 'q'
        options = sdpsettings('solver', 'quadprog', 'verbose', 1);
    elseif lower(optimizer_type) == 'f'
        options = sdpsettings('solver', 'fmincon', ...
                              'sdpa.maxIteration', 10, ...
                              'verbose', 0);                             
        %options = sdpsettings('solver', 'fmincon', ...
        %                    'relax', 0, ...
        %                    'warning', 0, ...
        %                    'showprogress', 1 ...
        %                    );
    elseif lower(optimizer_type) == 'o'
        options = sdpsettings('solver', 'OSQP', ... % Use OSQP solver
                      'verbose', 1, ...             
                      'osqp.max_iter', 30000, ...   % Set maximum iterations
                      'osqp.eps_abs', 1e-5, ...     % Absolute tolerance
                      'warmstart', 0);             % Disable warm start
    elseif lower(optimizer_type) == 'b'
        options = sdpsettings('solver', 'bmibnb', 'verbose', 1);
    else
        error('Error assigning solver options!');
    end
    
    diagnostics = optimize(constraints, objective, options);

    if diagnostics.problem == 0
                % g = value(g);
                u = value(u);
                x = value(x);
                y = value(y);
        if verbose
            disp("REFERENCE TRAJECTORY: "); disp(target');
            disp("DELTA: "); disp(value(delta)');
            disp("PEANALTY TERM: "); disp(value(objective));
        end
    else
        disp('Optimization failed:');
        disp(diagnostics.info);
        error('The problem did not solve successfully. YALMIP error code: %d', diagnostics.problem);
    end
end