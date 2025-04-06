function [A, B, C, D] = linearize(Fx, Gx, x, u, x_e, u_e)
    % Initialize symbolic Jacobians
    A_sym = [];
    B_sym = [];
    for i = 1:length(Fx)
        % Compute Jacobians for each dynamic equation in Fx
        A_sym = [A_sym; jacobian(Fx{i}, x)]; % State Jacobian
        B_sym = [B_sym; jacobian(Fx{i}, u)]; % Input Jacobian
    end

    % Measurement Jacobians
    C_sym = [];
    D_sym = [];
    for i = 1:length(Gx)
        % Compute Jacobians for each measurement equation in Gx
        C_sym = [C_sym; jacobian(Gx{i}, x)]; % State Jacobian
        D_sym = [D_sym; jacobian(Gx{i}, u)]; % Input Jacobian
    end

    % Ensure dimensions of substitution match
    if numel([x; u]) ~= numel([x_e; u_e])
        error('Mismatch in dimensions of variables and equilibrium values.');
    end

    % Evaluate Jacobians numerically at equilibrium
    A = double(subs(A_sym, [x; u], [x_e; u_e]));
    B = double(subs(B_sym, [x; u], [x_e; u_e]));
    C = double(subs(C_sym, [x; u], [x_e; u_e]));
    D = double(subs(D_sym, [x; u], [x_e; u_e]));
end