function output_dir = prepareOutputDir(mode)
    switch mode
        case 'plots'
            dirname = 'plots';
        case 'data'
            dirname = 'data';
        otherwise
            error('Unknown mode: %s', mode);
    end
    output_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'outputs', dirname);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
end
