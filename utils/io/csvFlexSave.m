function filename = csvFlexSave(prefix, varargin)
    %% Saves multiple sequences to a structured CSV file
    % Inputs:
    %   prefix: A string used as the base filename.
    %   varargin: A variable number of data inputs (can be numeric vectors or cells).
    
    if nargin < 2
        error('At least one data vector and a prefix must be provided.');
    end
    
    % Validate prefix
    if isempty(prefix) || ~ischar(prefix)
        error('The prefix must be a non-empty string.');
    end
    
    % Initialize structured data storage
    nInputs = numel(varargin);    
    d_csv = [];

    % Preprocess and validate each input
    for i = 1:nInputs
        col_i = varargin{i};
        col_i = col_i(:);    
        d_csv = [d_csv, col_i];
    end    

    % Construct the full filename (add .csv extension if not included)
    if ~endsWith(prefix, '.csv')
        filename = sprintf('%s.csv', prefix);
    else
        filename = prefix;
    end

    % Write data to CSV    
    try
        writematrix(d_csv, filename, 'Delimiter', ',', 'WriteMode', 'overwrite');
        fprintf('File saved to: <PATH TO MATLAB>\DDSF\outputs\data\%s', filename);
    catch matrixME
        fprintf('WRITEMATRIFailed to write data to CSV: %s', matrixME.message);
        try
            writecell(d_csv, filename, 'Delimiter', ',', 'WriteMode', 'overwrite');
        catch cellME
             error('Failed to write data to CSV: %s', cellME.message);
        end
    end
end