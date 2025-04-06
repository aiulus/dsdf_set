function sys = populate_system_struct(A, B, C, D, params)
    % Constructs the system struct with matrices, dimensions, and constraints.
    sys.A = A;
    sys.B = B;
    sys.C = C;
    sys.D = D;

    % Parameters and target
    sys.params = params;
    sys.target = params.target;
    sys.dims = getDims(A, B, C, D);
end

