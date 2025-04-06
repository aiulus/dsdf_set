function [u_d, y_d, x_d, u, y] = genDataNonlinear(lookup)
    %% Extract parameters
    sys = lookup.sys;
    dims = sys.dims;
    dt = sys.dt;

    Fh = sys.functions.Fh;
    Gh = sys.functions.Gh;

    params = lookup.sys.params;
    T = lookup.config.T;

    % State and input dimensions
    n = dims.n;
    m = dims.m;
    p = dims.p;

    %% Generate a random control input
    PE_input = inputSignalGenerator(lookup, T);

    % Initialize input-output storage
    u_d = zeros(m, T);
    y_d = zeros(p, T);
    x_d = zeros(n, T + 1);

    % Set initial state
    x_d(:, 1) = params.x_ini;

    %% Simulate the system
    for t = 1:T
        u_d(:, t) = PE_input(:, t); % Input at time t
    
        % Evaluate dynamics
        x_d(:, t + 1) = evalDynamicsFct(Fh, x_d(:, t), u_d(:, t), dt);
        y_d(:, t) = evalMsrmtFct(Gh,  x_d(:, t), u_d(:, t));
    end

    % Flatten the control inputs and outputs
    u = reshape(u_d, [], 1); % Reshapes into (T * m) x 1
    y = reshape(y_d, [], 1); % Reshapes into (T * p) x 1
end
