function [u, y] = getNearestEquilibrium(S_f, y_current, u_current)
    u_e = S_f.symbolic_solution.u_e;
    y_e = S_f.symbolic_solution.y_e;
    y = vectorFilter(y_current, y_e);
    u = vectorFilter(u_current, u_e);
end