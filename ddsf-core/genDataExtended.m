function data = genDataExtended(lookup)
    %% === Extract system and configuration ===
    sys = lookup.sys;
    A = sys.A; B = sys.B; C = sys.C; D = sys.D;

    m = lookup.dims.m;
    n = lookup.dims.n;
    p = lookup.dims.p;
    T = lookup.config.T;
    T_d = lookup.T_d;

    x0_mode = "single"; % default
    if isfield(lookup.data_options, "x0_mode")
        x0_mode = lookup.uncertainty_options.x0_mode;
    end

    %% === Generate PE input ===
    PE_input = inputSignalGenerator(lookup, T);
    u_d = zeros(m, T);

    if isfield(sys, "sets") && isfield(sys.sets, "X0") && isa(sys.sets.X0, 'zonotope')
        X0 = sys.sets.X0;
    else
        X0 = zonotope(sys.params.x_ini); % Fallback
    end

    %% === Sample-based mode ===
    if strcmp(x0_mode, "single")
        x0 = sys.params.x_ini;
            x_d = zeros(n, T+1);
            y_d = zeros(p, T);
            if lookup.opt_params.init
                x_d(:, 1) = x0;
                u_d(:, 1) = PE_input(:, 1);
                y_d(:, 1) = C * x_d(:, 1) + D * u_d(:, 1);
                low = 2;
            else
                x_d(:, 1) = x0;
                low = 1;
            end
            for i = low:T
                u_d(:, i) = PE_input(:, max(i - T_d, 1));
                x_d(:, i+1) = A * x_d(:, i) + B * u_d(:, i);
                y_d(:, i) = C * x_d(:, i) + D * u_d(:, i);
            end
            data = struct('mode', 'single', 'u_d', u_d, 'y_d', y_d, ...
                          'x_d', x_d, 'u', reshape(u_d, [], 1), 'y', reshape(y_d, [], 1));

    elseif strcmp(x0_mode, "sample")
        x0 = sample(X0);
        x_d = zeros(n, T+1);
        y_d = zeros(p, T);

        % Initialization
        if lookup.opt_params.init
            x_d(:, 1) = x0;
            u_d(:, 1) = PE_input(:, 1);
            y_d(:, 1) = C * x_d(:, 1) + D * u_d(:, 1);
            low = 2;
        else
            x_d(:, 1) = x0;
            low = 1;
        end

        for i = low:T
            if (i - T_d) < 1
                u_d(:, i) = 0;
            else
                u_d(:, i) = PE_input(:, i - T_d);
            end
            x_d(:, i+1) = A * x_d(:, i) + B * u_d(:, i);
            y_d(:, i) = C * x_d(:, i) + D * u_d(:, i);
        end

        data = struct( ...
            'mode', 'sample', ...
            'u_d', u_d, ...
            'y_d', y_d, ...
            'x_d', x_d, ...
            'u', reshape(u_d, [], 1), ...
            'y', reshape(y_d, [], 1) ...
        );

    %% === Zonotope propagation mode ===
    elseif strcmp(x0_mode, "set")
        X_zono = cell(1, T+1);
        Y_zono = cell(1, T);
        x_nominal = zeros(n, T+1);
        y_nominal = zeros(p, T);
        x_d = zeros(n, T+1);
        y_d = zeros(p, T);

        X_zono{1} = X0;
        x_nominal(:, 1) = getCenter(X0);

        for i = 1:T
            if (i - T_d) < 1
                u_d(:, i) = 0;
            else
                u_d(:, i) = PE_input(:, i - T_d);
            end
            X_zono{i+1} = A * X_zono{i} + B * u_d(:, i);
            Y_zono{i}   = C * X_zono{i} + D * u_d(:, i);

            x_nominal(:, i+1) = getCenter(X_zono{i+1});
            y_nominal(:, i)   = getCenter(Y_zono{i});

            % Also store one nominal trajectory for raw Hankel if needed
            x_d(:, i+1) = x_nominal(:, i+1);
            y_d(:, i)   = y_nominal(:, i);
        end

        data = struct( ...
            'mode', 'propagate', ...
            'u_d', u_d, ...
            'y_d', y_d, ...
            'x_d', x_d, ...
            'X_zono', {X_zono}, ...
            'Y_zono', {Y_zono}, ...
            'x_nominal', x_nominal, ...
            'y_nominal', y_nominal, ...
            'u', reshape(u_d, [], 1), ...
            'y', reshape(y_d, [], 1) ...
        );

    else
        error("Unknown data generation mode: '%s'. Use 'sample' or 'propagate'.", x0_mode);
    end
end
