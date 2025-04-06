%% Retrieves individual sequences from a CSV file
% Inputs:
%   filename: The path to the CSV file respective to the MATLAB directory.
%
% Outputs:
%   data: A cell array where each cell contains one sequence.
%   sequenceNames: A cell array of sequence names for reference.
%
% This function reads a CSV file saved with `csvFlexSave` and extracts individual
% sequences as numeric arrays.

function [data, configs] = csvFlexRead(filename)
    if nargin < 1 || isempty(filename) || ~ischar(filename)
        error('Filename must be a non-empty string.');
    end

    if ~isfile(filename)
        error('The file %s does not exist.', filename);
    end

    try
        fullPath = which(filename);
        if isempty(fullPath)
            % If `which` fails, try resolving relative path
            fullPath = fullfile(pwd, filename);
        end
        if ~isfile(fullPath)
            error('File not found. Ensure the file path is correct: %s', fullPath);
        end
    catch
        error('Error resolving full path for the file. Verify the input filename.');
    end

    % Read the CSV file into a table
    try
        data = readtable(filename, 'Delimiter', ',');
    catch err
        error('Error reading the CSV file: %s', err.message);
    end

    %data = data(:, 1:end-1);
    configs = table2cell(data(:, end));
    data = data(:, 1:end-1);
end
