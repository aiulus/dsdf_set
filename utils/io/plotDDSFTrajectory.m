function plotDDSFTrajectory(data)
    T = size(data.u_d, 2);
    t = 1:T;

    if strcmp(data.mode, "set") && isfield(data, "Y_zono")
        hold on;
        for i = 1:T
            plotFilledZonotope(data.Y_zono{i}, i);
        end
    end

    if isfield(data, "y_nominal")
        plot(t, data.y_nominal', 'LineWidth', 2);
    else
        plot(t, data.y_d', '--', 'LineWidth', 1.5);
    end

    xlabel("Time Step");
    ylabel("Output");
    grid on;
    legendStrings = arrayfun(@(i) sprintf("y_{%d}", i), 1:size(data.y_d,1), 'UniformOutput', false);
    legend(legendStrings{:}, 'Location', 'best');
end

function plotFilledZonotope(Z, t)
    % Loop over each dimension of output
    p = size(Z.Z, 1); % Number of output dimensions

    for j = 1:p
        Z_proj = project(Z, j); % Single dimension projection
        interval_bounds = interval(Z_proj);  % 2 x 1
        y_lb = interval_bounds(1);
        y_ub = interval_bounds(2);

        fill([t t], [y_lb y_ub], [0.8 0.8 1], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
end

