function gridPlotDeePC(configname, sys, sorted)
    % Extract data
    u_hist = sorted.u; 
    y_hist = sorted.y;

    % Prepare output directory
    output_dir = prepareOutputDir('plots');

    % Define plotting parameters
    T_sim = size(y_hist, 3); % Time steps
    time = 0:T_sim-1; % Time vector
    p = sys.dims.p; % Number of outputs (y)
    m = sys.dims.m; % Number of inputs (u)
    
    % ---- Plot Outputs (y_hist) ----
    plots_per_figure = 3; % Maximum number of plots per figure
    num_figures_y = ceil(p / plots_per_figure); % Number of figures for outputs
    
    for fig = 1:num_figures_y
        % Create a new figure for outputs
        figure('Position', [100, 100, 800, 600]); 
        tiledlayout(plots_per_figure, 1); 

        for sub = 1:plots_per_figure
            i = (fig - 1) * plots_per_figure + sub; % Global plot index
            if i > p, break; end % Stop if no more outputs to plot

            nexttile; hold on;

            % Plot y_hist
            plot(time, y_hist(i, :), 'b', 'LineWidth', 1.25, 'DisplayName', sprintf('y[%d]', i));

            % Add target line if available
            addTargetLine(sys.params.target, T_sim);

            % Configure plot
            title(sprintf('Output %d', i));
            xlabel('t');
            ylabel(sprintf('y[%d]', i));
            grid on;
            legend show;
            hold off;
        end

        % Save the figure
        saveAndClose(output_dir, sprintf('%s_outputs_fig%d', configname, fig));
    end

    % ---- Plot Inputs (u_hist) ----
    plots_per_figure = 3; % Maximum number of plots per figure
    num_figures_u = ceil(m / plots_per_figure); % Number of figures for inputs

    for fig = 1:num_figures_u
        % Create a new figure for inputs
        figure('Position', [100, 100, 800, 600]); 
        tiledlayout(plots_per_figure, 1); 

        for sub = 1:plots_per_figure
            i = (fig - 1) * plots_per_figure + sub; % Global plot index
            if i > m, break; end % Stop if no more inputs to plot

            nexttile; hold on;

            % Plot u_hist
            stairs(time, u_hist(i, :), 'b', 'LineWidth', 1.25, 'DisplayName', sprintf('u[%d]', i));

            % Add bounds if defined
            if isfield(sys, 'constraints') && isfield(sys.constraints, 'U')
                addBounds(time, sys.constraints.U(i, :), configname);
            end

            % Configure plot
            title(sprintf('Input %d', i));
            xlabel('t');
            ylabel(sprintf('u[%d]', i));
            grid on;
            legend show;
            hold off;
        end

        % Save the figure
        saveAndClose(output_dir, sprintf('%s_inputs_fig%d', configname, fig));
    end
end
