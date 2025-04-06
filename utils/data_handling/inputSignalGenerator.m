function u = inputSignalGenerator(lookup, length)
    % Generates input signals within (relaxed) system constraints.
    %
    % Parameters:
    %   sys               - System structure containing constraints and dimensions.
    %   mode              - String specifying the input signal type:
    %                       'rbs'            : Random Binary Signal
    %                       'scaled_rbs'     : Scaled Random Binary Signal
    %                       'scaled_gaussian': Scaled Gaussian Signal
    %   length             - Integer specifying the number of time steps.
    %   relaxation_factor - Scalar to relax constraints (e.g., 1.25 for 25% relaxation).
    %
    % Returns:
    %   u                 - Generated input signal matrix of size (m x length).

    sys = lookup.sys; 
    mode = lookup.data_options.datagen_mode;
    scale = lookup.data_options.scale;

    % Extract system dimensions and constraints
    m = sys.dims.m;
    lb = sys.constraints.U(:, 1);
    ub = sys.constraints.U(:, 2);

    % Replace infinite bounds with finite values
    lb(lb == -inf) = -1e+8;
    ub(ub == inf) = 1e+8;

    % Generate input signals based on the specified mode
    switch mode
        case 'prbs'
            % Random Binary Signal
            u = idinput([length, m], 'prbs', [0, 1], [-1, 1])';
        case 'sinusoid'
            u = customSinusoid(sys, length, scale, 1e2);
        case 'sinusoidal_sweep'
            u = sinusoidal_sweep(lookup, lb, ub, m, length);      
        case 'uniform'            
            u = uniform(lb, ub, m, length);
        case 'custom_uniform'
            u = customUniform(scale, lb, ub, m, length);
        case 'controlled_random'
            u = controlled_random(lookup, lb, ub, m, length);
        case 'white_noise'
            u = white_noise(lb, ub, m, length);
        otherwise
            error(['Unsupported mode: %s\n Please choose from ''prbs'', ''sinusoid'', ''sinusoidal_sweep'', ''uniform'', ' ...
                   '''custom_uniform'', ''controlled_random'', ''white_noise''.'], mode);
    end

    % Ensure signals are within relaxed bounds
    if lookup.data_options.safe
        u = max(min(u, ub), lb);
    end
end



