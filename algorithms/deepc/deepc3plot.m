function deepc3plot(time, y_sim, u_sim, sys)
    % Plot outputs
    figure(1);
    p = size(y_sim, 1);
    sys.target = sys.params.target;

    output_dir = prepareOutputDir('plots');

    for i=1:p
        subplot(p, 1, i);
        y_name = sprintf("y%d", i);
        hold on;
        plot(time, y_sim(1, :), 'r', 'LineWidth', 1.5, 'DisplayName', y_name); 
        plot(time, sys.constraints.Y(1)*ones(size(time)), 'b', 'LineStyle','-.', 'LineWidth', 1.5,'DisplayName', 'y_{min}')
        plot(time, sys.constraints.Y(2)*ones(size(time)), 'b', 'LineStyle','--', 'LineWidth', 1.5, 'DisplayName', '{y_{max}')
        xlabel('Iteration #');
        ylabel(sprintf('Output %d', i));

        target = sys.target * ones(1, size(time, 2));
        plot(time, target, 'm--', 'DisplayName', 'Target');
       
        grid on; legend show; hold off;
    end
    sgtitle("Control outputs over time");
    saveAndClose(output_dir, 'deepc-damper-outputs');
    
    % Plot inputs
    figure(2);
    m = size(u_sim, 1);

    for i=1:m
        subplot(m, 1, i);
        u_name = sprintf("u%d", i);
        hold on;
        stairs(time, u_sim(1, :), 'r', 'LineWidth', 1.5, 'DisplayName', u_name);
        xlabel('Iteration #');
        ylabel(sprintf('Control input %d', i));

        plot(time, sys.constraints.U(1)*ones(size(time)), 'b', 'LineStyle','-.', 'LineWidth', 1.5,'DisplayName', 'u_{min}')
        plot(time, sys.constraints.U(2)*ones(size(time)), 'b', 'LineStyle','--', 'LineWidth', 1.5,'DisplayName', 'u_{max}')

        grid on; legend show; hold off;
    end
    sgtitle("System inputs over time");
    saveAndClose(output_dir, 'deepc-damper-inputs');
end

