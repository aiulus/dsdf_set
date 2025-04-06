%% Equilibrium Computation
function S_f = getEquilibriumNonlinear(Fx, x, u)
    % Construct equilibrium equations
    eqns = [];
    for i = 1:length(Fx)
        eqns = [eqns; Fx{i} == 0];
    end
    
    % Solve for equilibrium
    vars = [x(:); u(:)]; % Combine state and input variables

    %% EXTRACT TRIVIAL SOLUTION
    trivial_sol = solve(eqns, vars, 'Real', true);

    % Check if the solution is empty
    if isempty(trivial_sol)
        error('No real equilibrium solution found.');
    end

    % Extract the solution into numeric arrays
    if isstruct(trivial_sol)
        % Single solution: Convert structure fields to array
        trivial_sol_vals = struct2array(trivial_sol); % Extract all variables in order
    elseif iscell(trivial_sol)
        % Multiple solutions: Use the first one
        trivial_sol_vals = struct2array(trivial_sol{1});
    else
        error('Unexpected solution format from solve.');
    end

    % Split into state and input equilibrium values
    num_states = numel(x); % Number of state variables
    num_inputs = numel(u); % Number of input variables
    trivial_x_e = double(trivial_sol_vals(1:num_states));   % Extract state equilibrium
    trivial_u_e = double(trivial_sol_vals(num_states+1:num_states+num_inputs)); % Extract input equilibrium

    % Reshape to match the dimensions of symbolic variables
    trivial_x_e = reshape(trivial_x_e, size(x));
    trivial_u_e = reshape(trivial_u_e, size(u));

    %% EXTRACT PARAMETRIZED SOLUTION
    
    sol = solve(eqns, vars, 'Real', false, 'ReturnConditions', true);

    % Check if the solution is empty
    if isempty(sol)
        error('No real equilibrium solution found.');
    end

    % Extract the solution into numeric arrays
    if isstruct(sol)
        % Single solution: Convert structure fields to array
        sol_vals = struct2array(sol); % Extract all variables in order
    elseif iscell(sol)
        % Multiple solutions: Use the first one
        sol_vals = struct2array({1});
    else
        error('Unexpected solution format from solve.');
    end

    x_e = sym(sol_vals(1:num_states));   % Extract state equilibrium
    u_e = sym(sol_vals(num_states+1:num_states+num_inputs)); % Extract input equilibrium

    % Reshape to match the dimensions of symbolic variables
    x_e = reshape(x_e, size(x));
    u_e = reshape(u_e, size(u));

    % Package equilibrium points into a struct
    S_f = struct( ...
            'trivial_solution', struct('x_e', trivial_x_e, 'u_e', trivial_u_e), ...
            'symbolic_solution', struct('x_e', x_e, 'u_e', u_e) ...
    );
end