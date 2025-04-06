function c_d = ddsfConservatismSingle(lookup, logs)
    sys = lookup.sys;
    lb = sys.constraints.Y(:, 1); 
    ub = sys.constraints.Y(:, 2);
    p = sys.dims.p;
    y_d = logs.y; yl_d = logs.yl;
    % T_sim = size(y_d, 2);
    c_i = zeros(1, p);
    for d=1:p    
        yi_d = y_d(d, :); yli_d = yl_d(d, :);
        ci_d = conservatism(yi_d, yli_d, lb(d), ub(d));
        c_i(d) = ci_d;
    end
    c_d = mean(c_i);
end