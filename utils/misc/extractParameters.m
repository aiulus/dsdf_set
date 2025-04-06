function params = extractParameters(mode, index, values, param_struct)
    % Extracts parameters for the current configuration based on the mode and index.

    switch mode
        case {'NvsTini', 'nt'}
            params = struct('T_ini', values(index, 1), 'N', values(index, 2), 'Q', param_struct.Q, 'R', param_struct.R);

        case {'QvsR', 'qr'}
            params = struct('Q', values(index, 1), 'R', values(index, 2), 'T_ini', param_struct.T_ini, 'N', param_struct.N);

        case 'mixed'
            [qr_idx, nt_idx] = ind2sub([size(param_struct.qr, 1), size(param_struct.nt, 1)], index);
            params = struct('T_ini', param_struct.nt(nt_idx, 1), 'N', param_struct.nt(nt_idx, 2), ...
                            'Q', param_struct.qr(qr_idx, 1), 'R', param_struct.qr(qr_idx, 2));
    end
end