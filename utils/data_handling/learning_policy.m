%   Can later be changed to:
%       u_l = learning_policy(y, y_d)
%
%   INPUTS:
%       y   - Current system output (px1)
%       y_d - Desired system output    

function u_l = learning_policy(lookup)
    sys = lookup.sys;
    m = lookup.sys.dims.m;
    % L = lookup.config.N + 2 * lookup.config.T_ini;
    Np = lookup.config.N;
    mode = lookup.data_options.datagen_mode;
    scale = lookup.data_options.scale;

    U = sys.constraints.U; 
    U_relaxed = updateBounds(U, scale);    
    
    lb = U_relaxed(:, 1);
    lb(lb == -inf) = -1e+8;
    ub = U_relaxed(:, 2);
    ub(ub == inf) = 1e+8;
    
    switch mode
        case 'prbs'
            u_l = idinput([m, Np], 'prbs', [0, 1], [-1,1]);
        case 'sinusoid'
            u_l = customSinusoid(lookup.sys, lookup.T_sim, scale, 1e+8);
         case 'sinusoidal_sweep'
            u_l = sinusoidal_sweep(lookup, lb, ub, m, Np);       
        case 'uniform'
            u_l = uniform(lb, ub, m, Np);
        case 'custom_uniform'
            u_l = customUniform(scale, lb, ub, m, Np);
        case 'controlled_random'
            u_l = controlled_random(lookup, lb, ub, m, Np);
        case 'white_noise'
            u_l = white_noise(lb, ub, m, Np);
    end
end




