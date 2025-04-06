function [u_d, y_d, x_d, u, y] = simulateNonlinearWithDelay(sys, config)
    %% Extract Parameters
    params = sys.params;
    Fh = sys.functions.Fh; % Nonlinear dynamics handle
    Gh = sys.functions.Gh; % Measurement function handle
    dt = params.dt;        % Sampling time
    T = config.T;          % Simulation time horizon
    T_d = ceil(params.time_delay / dt); % Time delay in steps

    % Dimensions
    n = sys.dims.n; % State dimension
    m = sys.dims.m; % Input dimension
    p = sys.dims.p; % Output dimension

    %% Generate Random Input
    PE_input = inputSignalGenerator(sys, T); % Persistent excitation input

    % Initialize storage
    u_d = zeros(m, T);   % Delayed control input
    y_d = zeros(p, T);   % Output
    x_d = zeros(n, T+1); % State trajectory
    u = zeros(m * T, 1); % Flattened input
    y = zeros(p * T, 1); % Flattened output

    % Set initial state
    x_d(:, 1) = params.x_ini;

    %% Simulation Loop
    for t = 1:T
        % Handle delay: apply delayed input or zero if within delay period
        if t > T_d
            u_d(:, t) = PE_input(:, t - T_d); % Use delayed input
        else
            u_d(:, t) = zeros(m, 1); % No input during delay period
        end

        % Evaluate nonlinear dynamics
        x_d(:, t+1) = evalDynamicsFct(Fh, x_d(:, t), u_d(:, t), dt);

        % Evaluate measurement output
        y_d(:, t) = evalMsrmtFct(Gh, x_d(:, t), u_d(:, t));
    end

    % Flatten input and output for compatibility with DDSF frameworks
    u = reshape(u_d, [], 1); % Reshaped into (T * m) x 1
    y = reshape(y_d, [], 1); % Reshaped into (T * p) x 1
end
