function [Ad, Bd, Cd, Dd] = discretize_system(A, B, C, D, dt)
    % Discretizes a continuous-time state-space system.
    sys_cont = ss(A, B, C, D);
    sys_disc = c2d(sys_cont, dt);
    [Ad, Bd, Cd, Dd] = ssdata(sys_disc);
end