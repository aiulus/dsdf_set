function sys = constraint_handler(sys, params)
    dims = sys.dims;
    largeval = 1e+8;

    % Assign constraints
    if isfield(params, 'u_min')
        if max(size(params.u_min)) == 1
            u_min = repmat(params.u_min, dims.m, 1);
        else
            u_min = params.u_min;
        end
    else
        u_min = repmat(-largeval, dims.m, 1);
    end

    if isfield(params, 'u_max')
        if max(size(params.u_max)) == 1
            u_max = repmat(params.u_max, dims.m, 1);
        else
            u_max = params.u_max;
        end
    else
        u_max = repmat(largeval, dims.m, 1);
    end

    if isfield(params, 'y_min')
        if max(size(params.y_min)) == 1
            y_min = repmat(params.y_min, dims.p, 1);
        else
            y_min = params.y_min;
        end
    else
        y_min = repmat(-largeval, dims.p, 1);
    end

    if isfield(params, 'y_max')
        if max(size(params.y_max)) == 1
            y_max = repmat(params.y_max, dims.p, 1);
        else
            y_max = params.y_max;
        end
    else
        y_max = repmat(largeval, dims.p, 1);
    end

    sys.constraints.U = [u_min, u_max];
    sys.constraints.Y = [y_min, y_max];
end

