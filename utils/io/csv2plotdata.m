function [plotdata, configs] = csv2plotdata(filename, mode, T_sim, dims)
    m = dims.m; p = dims.p;
    raw_data = readmatrix(filename, 'Delimiter', ',');
    configs = raw_data(:, end); % Last column contains configurations
    raw_data = raw_data(:, 1:end-1); % All but the last column are numeric

    % Preallocate structures

    nruns = size(raw_data, 1);
    plotdata = struct('u', [], 'ul', [], 'y', [], 'yl', []);

    % Parse data based on mode
    switch mode
        case {'ddsf', 'u-ddsf'}
            % Extract u and ul for DDSF (dimensions m)
            for i = 1:nruns
                start_u = 1;
                end_u = T_sim * m;
                plotdata.u(i, :, :) = reshape(raw_data(i, start_u:end_u), [m, T_sim]);

                start_ul = end_u + 1;
                end_ul = start_ul + T_sim * m - 1;
                plotdata.ul(i, :, :) = reshape(raw_data(i, start_ul:end_ul), [m, T_sim]);
            end

        case 'y-ddsf'
            % Extract y and yl for DDSF (dimensions p)
            for i = 1:nruns
                start_y = 1;
                end_y = T_sim * p;
                plotdata.y(i, :, :) = reshape(raw_data(i, start_y:end_y), [p, T_sim]);

                start_yl = end_y + 1;
                end_yl = start_yl + T_sim * p - 1;
                plotdata.yl(i, :, :) = reshape(raw_data(i, start_yl:end_yl), [p, T_sim]);
            end

        case 'deepc'
            % Extract u (m) and y (p) for DeePC
            for i = 1:nruns
                start_u = 1;
                end_u = T_sim * m;
                plotdata.u(i, :, :) = reshape(raw_data(i, start_u:end_u), [m, T_sim]);

                start_y = end_u + 1;
                end_y = start_y + T_sim * p - 1;
                plotdata.y(i, :, :) = reshape(raw_data(i, start_y:end_y), [p, T_sim]);
            end

        otherwise
            error('Unsupported mode: %s', mode);
    end
end