function save_plots(output_dir, suffix)
    figures = findall(0, 'Type', 'figure');
    for i = 1:numel(figures)
        fig_name_base = sprintf('plot_%s_fig%d', suffix, i);
        fig_path = fullfile(output_dir, [fig_name_base, '.tex']);

        % Increment file name if it exists
        counter = 1;
        while exist(fig_path, 'file')
            fig_name_base = sprintf('plot_%s_fig%d_v%d', suffix, i, counter);
            fig_path = fullfile(output_dir, [fig_name_base, '.tex']);
            counter = counter + 1;
        end

        % Save the figure as .tex
        matlab2tikz(fig_path, 'figurehandle', figures(i));
        fprintf('Saved figure: %s\n', fig_path);
    end
end
