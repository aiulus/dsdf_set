function sys = populate_system_struct(A, B, C, D, params)
    % Constructs the system struct with matrices, dimensions, and constraints.
    sys.Ac = A;
    sys.Bc = B;
    sys.Cc = C;
    sys.Dc = D;

    % Dimensions
    sys.dims.state = size(A, 1);
    sys.dims.input = size(B, 2);
    sys.dims.output = size(C, 1);

    % Constraints
    sys.constraints.U = [params.u_min, params.u_max];
    sys.constraints.Y = [params.y_min, params.y_max];

    % Parameters and target
    sys.params = params;
    sys.target = params.target;

    [Ad, Bd, Cd, Dd] = discretize_system(A, B, C, D, params.dt);

    %% Store the discretized system as (A, B, C, D)
    sys.A = Ad;
    sys.B = Bd;
    sys.C = Cd;
    sys.D = Dd;
end