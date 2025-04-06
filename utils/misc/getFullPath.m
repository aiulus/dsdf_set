function fullpath = getFullPath(prefix)    
    % Define the base output directory
    output_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'outputs', 'data');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % Ensure the directory exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir); % Create the directory if it doesn't exist
    end

    % Construct prefixes for filenames
    % prefix_u = fullfile(base_output_dir, prefix);
    fullpath = fullfile(output_dir, prefix);
end