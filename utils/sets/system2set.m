function sys = system2set(sys)
    try
        n = size(sys.A, 1);
        m = size(sys.B, 2);
        p = size(sys.C, 1);
    
        % Initial condition set
        x0 = sys.params.x_ini;
        if isempty(x0) || any(isnan(x0)) || any(isinf(x0))
            x0 = zeros(n, 1); % default
        end
        x0 = double(x0);
        if isrow(x0), x0 = x0'; end
    
        %% TODO: add epsilon to global parameters
        epsilon = 0.01;
        sys.sets.X0 = safeZonotope(x0, epsilon * eye(n)); % Small uncertainty ball

        %% Input/Output Constraints
        umin = double(sys.params.u_min(:));
        umax = double(sys.params.u_max(:));
        ymin = double(sys.params.y_min(:));
        ymax = double(sys.params.y_max(:));
    
        sys.sets.constraints.U = zonotope(interval(umin, umax));
        sys.sets.constraints.Y = zonotope(interval(ymin, ymax));
    
        % Terminal safe set as zonotope
        if isfield(sys.S_f, "trivial_solution")
            x_eq = double(sys.S_f.trivial_solution.x_e(:));
            if isrow(x_eq), x_eq = x_eq'; end
            if isnumeric(x_eq)
                sys.sets.S_f_zono = safeZonotope(x_eq, 1e-3 * eye(n));
            end
        end
    catch ME
        fprintf('Zonotope conversion skipped: %s', ME.message);
    end
end