function x_next = evalDynamicsFct(Fh, x_current, u_current, dt)
    x_next = zeros(size(x_current));
    n = numel(Fh);
    for i=1:n
        fi = Fh{i};
        fi_val = fi(x_current, u_current);
        x_next_i = x_current(i) + dt * fi_val;
        x_next(i) = x_next_i;
    end
end