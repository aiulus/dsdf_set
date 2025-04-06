function u = customSinusoid(sys, T_sim, scale, inf_sub)
    % Generates a sinusoidal input signal with increasing amplitude within system constraints.
    %
    % Parameters:
    %   sys     - System structure containing constraints and dimensions.
    %   T_sim   - Integer specifying the number of time steps.
    %   scale   - Scalar to relax constraints (e.g., 1.25 for 25% relaxation).
    %   inf_sub - Substitute value for infinite boundaries (e.g., 1e+6).
    %
    % Returns:
    %   u       - Generated input signal matrix of size (m x T_sim).
    
    dt = sys.params.dt; % Time step [s]
    n_cycles = 5; % Define the desired number of cycles
    freq = n_cycles / (T_sim * dt); % Calculate frequency
    
    % Extract system dimensions and constraints
    m = sys.dims.m;
    lb = sys.constraints.U(:, 1);
    ub = sys.constraints.U(:, 2);
    
    % Replace infinite bounds with specified substitute value
    lb(lb == -inf) = -inf_sub;
    ub(ub == inf) = inf_sub;
    
    % Apply relaxation factor to constraints
    lb = lb * scale;
    ub = ub * scale;
    
    % Time vector
    t = (0:T_sim-1) * dt;
    
    % Initialize input signal matrix
    u = zeros(m, T_sim);
    
    % Generate sinusoidal signals with increasing amplitude
    for i = 1:m
        % Linearly increasing amplitude from a small value to the full range
        initial_amplitude = 0.1 * (ub(i) - lb(i)) / 2; % Start at 10% of half the range
        final_amplitude = (ub(i) - lb(i)) / 2;        % Max amplitude to stay within bounds
        amplitude = linspace(initial_amplitude, final_amplitude, T_sim);
        
        % Calculate midpoint for centering the sinusoid
        midpoint = (ub(i) + lb(i)) / 2;
        
        % Generate sinusoidal signal centered within bounds
        u(i, :) = midpoint + amplitude .* sin(2 * pi * freq * t);
        
        % Ensure the signal stays within the relaxed bounds
        u(i, :) = max(min(u(i, :), ub(i)), lb(i));
    end
end
