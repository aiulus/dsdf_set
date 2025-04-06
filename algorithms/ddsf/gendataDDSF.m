function [u_d, y_d, x_d, u, y] = gendataDDSF(lookup)
    %% Extract parameters
    sys = lookup.sys;
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;

    m = lookup.dims.m;
    n = lookup.dims.n;
    p = lookup.dims.p;

    T = lookup.config.T;

    %% Generate a random control input
    PE_input = inputSignalGenerator(lookup, T);

    % Initialize input-output storage
    u_d = zeros(m, T);
    y_d = zeros(p, T);
    x_d = zeros(n, T + 1);

    % Generate data by simulating the system on random inputs for L steps
    low = 1;
    if lookup.opt_params.init
        x_d(:, 1) = lookup.sys.params.x_ini;
        y_d(:, 1) = C * x_d(:, 1) + D * u_d(:, 1);
        u_d(:, 1) = PE_input(:, 1);
        low = 2;
    end
    
    T_d = lookup.T_d;
    for i = low:T
        if (i - T_d) < 1
            u_d(:, i) = 0;
        else
            u_d(:, i) = PE_input(:, i - T_d);
        end        
        x_d(:, i + 1) = A * x_d(:, i) + B * u_d(:, i);
        y_d(:, i) = C * x_d(:, i) + D * u_d(:, i);
    end

    % Flatten the control inputs and outputs
    u = reshape(u_d, [], 1); % Reshapes into (T * sys.params.m) x 1
    y = reshape(y_d, [], 1); % Reshapes into (T * sys.params.p) x 1

end

