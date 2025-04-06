%% Saves up to three individual sequences
function save2csv(col1, col2, col3, prefix)
    if nargin < 3
        error('At least, col1, and col2, and prefix must be provided.');
    end

    col1 = col1(:);
    col2 = col2(:);

    if isempty(prefix)
        prefix = '';
    end
    if nargin < 4 || isempty(col3)
        d_csv = [col1, col2];
    else
        col3 = col3(:);
        d_csv = [col1, col2, col3];
    end

    %output_dir = 'outputs/data';
    output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'outputs', 'data');

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    filename = fullfile(output_dir, sprintf('%s.csv', prefix));
    try
        writematrix(d_csv, filename, 'Delimiter', ',', 'WriteMode', 'overwrite');
    catch matrixME
        fprintf(matrixME.message);
        try
            writecell(d_csv, filename, 'Delimiter', ',', 'WriteMode', 'overwrite');
        catch cellME
            fprintf(cellME.message);
        end
    end
end