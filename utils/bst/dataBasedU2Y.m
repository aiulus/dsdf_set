%% TODO: Must use an initial trajectory of length T_ini, 
%%       Hankel matrices of order (T_ini + 1)
function y_pred = dataBasedU2Y(u_l, u_ini, y_ini, H_u, H_y)
    if iscell(H_u)
        H_u = cell2mat(H_u); % Convert to a numeric array
    end
    
    if iscell(H_y)
        H_y = cell2mat(H_y); % Convert to a numeric array
    end
    
    numcols = size(H_u, 2);    
    p = max(size(y_ini));

    if size(H_y, 2) ~= numcols
        error('The Hankel matrices H_u and H_y must have a matching number of columns!');
    end

    H = [H_u; H_y];

    alpha = sdpvar(numcols, 1);
    y_pred = sdpvar(p, 1);

    u = [u_ini; u_l];
    y = [y_ini; y_pred];
    traj = [u; y];

    constraints = traj == H * alpha;
    objective = alpha' * alpha;

    options = sdpsettings('solver', 'quadprog', 'verbose', 0);
    diagnostics = optimize(constraints, objective, options);
    
    if diagnostics.problem == 0
        % Extract the optimal prediction
        y_pred = value(y_pred);
    else
        error('Optimization problem could not be solved! Solver message: %s', diagnostics.info');
    end
   
end