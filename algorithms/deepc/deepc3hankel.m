function H = deepc3hankel(lookup)
    % Extract parameters
    n = lookup.dims.n;
    m = lookup.dims.m;
    p = lookup.dims.p;
    T = lookup.config.T;
    L = lookup.config.L;
    T_ini = lookup.config.T_ini;
    N = lookup.config.N;
    
    % Check dimensions
    if T <= L
        error('Insufficient data length: T must be greater than L.');
    end

    for i = 0:L-1
        U(m*i+1:m*(i+1), :) = lookup.data.u_data(:, i+1:T-L+i+1);
        Y(p*i+1:p*(i+1), :) = lookup.data.y_data(:, i+1:T-L+i+1);
        X(n*i+1:n*(i+1), :) = lookup.data.x_data(:, i+1:T-L+i+1);
    end
    
    full_rank = PEness_check(U);
    isPE = deepc_PEness_check(lookup.data.u_data, T_ini, N, lookup.sys);
    isPE = isPE & full_rank;

    if ~isPE
        error('Persistency of excitation check failed. Please provide richer input data or adjust T_ini and N.');
    end

    Up = U(1:m*T_ini, :);
    Uf = U(end -(m*N)+1:end, :);
    Yp = Y(1:p*T_ini, :);
    Yf = Y(end - (p*N)+1:end, :);

    H = struct( ...
        'U', U, ...
        'Y', Y, ...
        'X', X, ...
        'Up', Up, ...
        'Yp', Yp, ...
        'Uf', Uf, ...
        'Yf', Yf ...
        );
end

