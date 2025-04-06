function [x_data, y_data, u_data] = deepc3generateData(sys, dims, run_config)
    % rng(0, 'twister');
    % Initialize containers for data
    u_data = zeros(dims.m, run_config.T);
    y_data = zeros(dims.p, run_config.T);
    x_data = zeros(dims.n, run_config.T + 1);

    % x = rand(dims.n, 1); % Random initial state
    x = sys.params.x_ini;
    x_data(:, 1) = x;
    
    for t = 1:run_config.T
        u = rand(dims.m, 1);
        y = sys.C* x + sys.D * u;
        x = sys.A * x + sys.B * u;
        u_data(:, t) = u;
        y_data(:, t) = y;
        x_data(:, t + 1) = x;
    end
end