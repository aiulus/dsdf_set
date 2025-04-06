function [values, nruns, param_struct] = getModeConfig(mode, vals)
    % Maps the mode to the appropriate parameter values and configurations.

    switch mode
        case {'NvsTini', 'nt'}
            values = vals.NvsTini;
            nruns = size(values, 1);
            param_struct = struct('Q', -1, 'R', -1);

        case {'QvsR', 'qr'}
            values = vals.QvsR;
            nruns = size(values, 1);
            param_struct = struct('T_ini', -1, 'N', -1);

        case 'mixed'
            values = vals.mixed;
            nruns = size(values.qr, 1) * size(values.nt, 1);
            param_struct = struct('qr', values.qr, 'nt', values.nt);

        otherwise
            error('Unsupported mode "%s".', mode);
    end
end