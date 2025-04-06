systype = "cruise_control";
rng(0, 'twister'); % Set seed and generator
T_sim = 25;
%[mpc_t, mpc_u, mpc_y] = runMPC(systype, T_sim);
[mpc_t, mpc_u, mpc_y] = mpcMainF(systype, T_sim);
[deepc_t, deepc_u, deepc_y] = runDeePC(systype, T_sim);
plotMPCvsDeePC(systype, mpc_t, mpc_u, mpc_y, deepc_t, deepc_u, deepc_y);

function plotMPCvsDeePC(systype, mpc_t, mpc_u, mpc_y, deepc_t, deepc_u, deepc_y)
    sys = deepc2_systems(systype);
    % Plot outputs
    figure(1);
    p = size(mpc_y, 1);
    
    for i=1:p
        subplot(p, 1, i);
        y_name = sprintf("y%d", i);
        hold on;
        plot(mpc_t, mpc_y, 'r', 'LineWidth', 1.25, 'DisplayName', strcat("MPC - ", y_name)); 
        plot(deepc_t, deepc_y, 'b', 'LineWidth', 1.25, 'LineStyle', ':', 'DisplayName', strcat("DeePC - ", y_name)); 
        xlabel('Iteration #');
        ylabel(sprintf('Output %d', i));

        target = sys.params.target * ones(1, size(mpc_t, 2));
        plot(mpc_t, target, 'm--', 'DisplayName', 'Target');
       
        grid on; legend show; hold off;
    end
    sgtitle("MPC vs. DeePC: Control inputs over time");
    
    % Plot inputs
    figure(2);
    m = size(mpc_u, 1);

    for i=1:m
        subplot(m, 1, i);
        u_name = sprintf("u%d", i);
        hold on;
        stairs(mpc_t, mpc_u, 'r', 'LineWidth', 1.25, 'DisplayName', strcat("MPC - ", u_name));
        stairs(deepc_t, deepc_u, 'b', 'LineWidth', 1.25, 'LineStyle', ':', 'DisplayName', strcat("DeePC - ", u_name));
        xlabel('Iteration #');
        ylabel(sprintf('Control input %d', i));

        % Plot boundaries
        if bounds(1) ~= -inf
            plot(mpc_t, bounds(1) * ones(size(mpc_t)), 'm--', 'DisplayName', 'Lower Bound');
        end
        if bounds(2) ~= inf
            plot(mpc_t, bounds(2) * ones(size(mpc_t)), 'k--', 'DisplayName', 'Upper Bound');
        end
        grid on; legend show; hold off;
    end
    sgtitle("MPC vs. DeePC: Control inputs over time");
    figure(1);
    matlab2tikz('outputs/plots/new-cc-deepc_vs_mpc-outputs.tex');
    figure(2);
    matlab2tikz('outputs/plots/new-cc-deepc_vs_-inputs.tex');
end

