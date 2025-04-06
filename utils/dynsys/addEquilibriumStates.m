%% TODO: NOT IN USE !!! - REMOVE 
function sys = addEquilibriumStates(sys)
    U = sys.constraints.U; % Input constraints set
    Y = sys.constraints.Y; % Output constraints set

    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;

    n = sys.params.n;
    k = min(size(B)); % Number of SVD values needed for the pseudo-inverse

    % x_s = sdpvar(n, 1);

    T_u = sparse_pinv(B, k) * (eye(n) - A);
    T_y = C + D * T_u;
    T_y_tilde = C*sparse_pinv((eye(n) - A), rank((eye(n) - A)))*B - D;
    
    % Compute the column space of T_u
    %[U_eq, ~, ~] = svd(T_u, 'econ');
    %[Y_eq, ~, ~] = svd(T_y, 'econ');

    % Represent the safe sets with polytopes
    U_s = populateSafeSet(T_u, U);
    Y_s = populateSafeSet(T_y, Y);

    sys.S_f = struct( ...
        'U_s', U_s, ...
        'Y_s', Y_s ...
        );
end



function X_s = populateSafeSet(T_x, X)
    numDataPoints = 1e5;
    dim = size(T_x, 2);
    lb = X(:, 1);
    ub =  X(:, 2); 
    X_s = [];
    t = 1;

    while t < numDataPoints + 1
        % Random Gaussian vector with increased variance, multiplied by a
        % random magnitude
        z = randn(1,1) * idinput([dim, 1], 'rgs', [0, 1], [-1,10]).';
        x = T_x * z';
        within_range = all(x >= lb & x <= ub);
        
        if within_range
            X_s(:, t) = x;
            t = t + 1;
        end
    end
end




