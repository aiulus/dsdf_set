function saveAndClose(output_dir, configname)
    saveas(gcf, fullfile(output_dir, strcat(configname, '_plot.png')));
    try 
        cleanfigure();
        matlab2tikz(fullfile(output_dir, strcat(configname, '_plot.tex')), 'standalone', true); 
    catch
        warning('Error exporting to TikZ'); 
    end
    close(gcf);
end