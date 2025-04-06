function S_f = setEquilibriaDDSF(sys)
    A = sys.A; B = sys.B; C = sys.C; D = sys.D;

    % Define dimensions
    n = size(A, 1); % Number of states
    m = size(B, 2); % Number of inputs

    % Dynamically create symbolic vectors for states and inputs
    x = sym('x', [n, 1]); % State vector of size n x 1
    u = sym('u', [m, 1]); % Input vector of size m x 1

    % Define system dynamics
    dx = A*x + B*u;

    % Solve for equilibrium
    trivial_sol = struct2array(solve(dx == 0, [x; u], 'Real', true));
    x_e_triv = trivial_sol(:, 1:n).';
    u_e_triv = trivial_sol(:, n+1:end).';

    % Compute output at equilibrium
    y_e_triv = C*x_e_triv + D*u_e_triv;

    sol = struct2array(solve(dx == 0, [x; u], 'Real', false, 'ReturnConditions', true));
    x_e = sym(sol(1:n).');
    u_e = sym(sol(n+1:n+m).');
    y_e = C*x_e + D*u_e;

    % Package equilibrium points into a struct
    S_f = struct( ...
            'trivial_solution', struct('x_e', x_e_triv, 'u_e', u_e_triv, 'y_e', y_e_triv), ...
            'symbolic_solution', struct('x_e', x_e, 'u_e', u_e, 'y_e', y_e) ...
    );
end


