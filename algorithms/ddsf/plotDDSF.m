function plotDDSF(time, logs, lookup)    
    sys = lookup.sys;
    ul_hist = logs.ul_t; % s
    u_hist = logs.u;
    y_hist = logs.y;
    yl_hist = logs.yl;
    output_dir = prepareOutputDir('plots');

    % --- Plot Learning Inputs vs Safe Inputs ---        
    figure(1);
    m = sys.dims.m;
    plots_per_figure = 3; % Max number of plots per figure
    num_figures = ceil(m / plots_per_figure); % Number of figures needed
    
    for fig = 1:num_figures
        figure('Position', [100, 100, 800, 600]); % Fixed figure size
                tiledlayout(plots_per_figure, 1); % Max 3 subplots per figure
    
            for sub = 1:plots_per_figure
                i = (fig - 1) * plots_per_figure + sub; % Global index for subplot
                if i > m, break; end % Stop if no more subplots are needed

                nexttile; grid on;
                stairs(time, ul_hist(i, :), 'r', 'LineStyle', ':','LineWidth', 1.75, 'DisplayName', sprintf('ul[%d]', i));
                hold on;
                stairs(time, u_hist(i, :), 'b', 'LineWidth', 1.25, 'DisplayName', sprintf('u[%d]', i));
            
                bounds = sys.constraints.U(i, :);
            
                % Plot boundaries
                if bounds(1) ~= -inf
                    plot(time, bounds(1) * ones(size(time)), 'b', 'LineStyle', '-.', 'DisplayName', 'Lower Bound');
                end
                if bounds(2) ~= inf
                    plot(time, bounds(2) * ones(size(time)), 'b', 'LineStyle', '--','DisplayName', 'Upper Bound');
                end
            
                title(sprintf('Learning vs Safe Input %d', i));
                xlabel('t');
                ylabel(sprintf('Input %d', i));
                grid on;
                legend show;
                hold off;
            end
            sgtitle('Learning Inputs vs. Safe Inputs');

             % Save the current figure
            prefix = sprintf('U-ddsf-%s-fig%d', lookup.systype, fig);
            saveas(gcf, fullfile(output_dir, strcat(prefix, '-plot.png')));
            matlab2tikz(fullfile(output_dir, strcat(prefix, '.tex')));
            close(gcf); 
    end        

    % --- Plot Resulting Outputs ---    
    figure(2);
    p = sys.dims.p;
    plots_per_figure = 3; % Max number of plots per figure
    num_figures = ceil(p / plots_per_figure);

    for fig = 1:num_figures
        figure('Position', [100, 100, 800, 600]); % Set figure size
        tiledlayout(plots_per_figure, 1); % Max plots per figure

        for sub = 1:plots_per_figure
            i = (fig - 1) * plots_per_figure + sub; % Global plot index
            if i > p, break; end % Stop if we exceed the total number of plots
    
            nexttile; hold on;
    
            % Evaluate dynamic bounds            
            factor = 3; delta = (factor - 1) / 2;
            lb = bounds(1); ub = bounds(2);        
            lower_bound = lb - delta .* abs(lb);
            upper_bound = ub + delta .* abs(ub);

            % Filter out-of-bound values
            yl_temp = yl_hist(i, :); % Temporary copy for filtering
            out_of_bounds = (yl_temp < lower_bound) | (yl_temp > upper_bound); % Out-of-bound indices
            yl_temp(out_of_bounds) = NaN; % Replace out-of-bound values with NaN for plotting
            
            % Plot `y_hist`
            stairs(time, y_hist(i, :), 'b', 'LineWidth', 1.75, 'DisplayName', sprintf('y[%d]', i));
            
            % Plot `yl_hist` with filtered values
            stairs(time, yl_temp, 'r', 'LineStyle', ':', 'LineWidth', 1.25, 'DisplayName', sprintf('yl[%d]', i));
    
            % Add markers for out-of-bound points
            if any(out_of_bounds)
                % Indices of out-of-bound points
                out_of_bounds_idx = find(out_of_bounds);
    
                % Mark the lower bound breaches
                below_bounds_idx = out_of_bounds_idx(yl_hist(i, out_of_bounds_idx) < lower_bound);
                if ~isempty(below_bounds_idx)
                    plot(time(below_bounds_idx), lower_bound * ones(size(below_bounds_idx)), ...
                        'kv', 'MarkerFaceColor', 'r', 'DisplayName', 'Out of Bounds (Below)');
                end
    
                % Mark the upper bound breaches
                above_bounds_idx = out_of_bounds_idx(yl_hist(i, out_of_bounds_idx) > upper_bound);
                if ~isempty(above_bounds_idx)
                    plot(time(above_bounds_idx), upper_bound * ones(size(above_bounds_idx)), ...
                        'k^', 'MarkerFaceColor', 'g', 'DisplayName', 'Out of Bounds (Above)');
                end
            end
    
            if bounds(1) ~= -inf
                plot(time, bounds(1) * ones(size(time)), 'm--', 'DisplayName', 'Lower Bound');
            end
            if bounds(2) ~= inf
                plot(time, bounds(2) * ones(size(time)), 'm--', 'DisplayName', 'Upper Bound');
            end
    
            % Annotate dynamic bounds and optional targets
            if lookup.opt_params.target_penalty
                pi = sys.params.target(i);
                plot(time, pi * ones(size(time)), 'g--', 'DisplayName', 'Target');
            end
        end
        sgtitle("Resulting Outputs");

        % Save the current figure
        prefix = sprintf('Y-ddsf-%s-fig%d', lookup.systype, fig);
        saveas(gcf, fullfile(output_dir, strcat(prefix, '-plot.png')));
        matlab2tikz(fullfile(output_dir, strcat(prefix, '.tex')));
        close(gcf); 
    end
end
