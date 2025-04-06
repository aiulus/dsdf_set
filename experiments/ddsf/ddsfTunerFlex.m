function [u, ul, y, yl, descriptions, filename] = ddsfTunerFlex(mode, vals, systype, T_sim, toggle_save)
    % DDSF Tuner: A function to optimize parameters for DDSF.
    % 
    % INPUTS:
    %   mode        - Tuning mode ('r', 'NvsTini', 'constraints', 'mixed').
    %   vals        - Struct containing tuning parameter configurations.
    %   systype     - System type for DDSF.
    %   T_sim       - Simulation time.
    %   toggle_save - Boolean flag to save results (default: false).
    %
    % OUTPUTS:
    %   u, ul       - Resultant control inputs.
    %   y, yl       - Resultant system outputs.
    %   descriptions - Descriptions of each tuning configuration.

    % Default value for toggle_save
    if nargin < 5, toggle_save = true; end

    % Validate inputs
    validateInputs(mode, vals);

    max_tries = 5; % Allows for re-runs under the current configuration

    % Resolve mode and parameter configurations
    [values, nruns, param_settings] = resolveMode(mode, vals, max_tries);

    % Initialize output variables
    u = cell(1, nruns); ul = cell(1, nruns); 
    y = cell(1, nruns); yl = cell(1, nruns);
    descriptions = cell(1, nruns);

    % Create output directory if needed
    output_dir = fullfile('..', 'outputs', 'plots', 'ddsf_tuner');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % Main loop for configurations
    tic;
    for i = 1:nruns
        [param_desc, param_values] = buildDescription(mode, param_settings, values, i, systype, T_sim);
        fprintf('("------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);

        % Display progress for elapsed time
        if i > 1, timeEstimator(toc, i, nruns); end

        % Attempt to run DDSF
        for k = 1:param_settings.max_tries
            try
                [u{i}, ul{i}, y{i}, yl{i}, descriptions{i}] = executeDDSF(systype, T_sim, param_values, param_desc);
                break;
            catch ME
                fprintf(['Attempt to run runDDSF.m (conf.: %s) failed at: %s\n ' ...
                    'Message: = %s\n. Trying again...\n'], param_desc, ME.stack(1).name, ME.message);
            end
        end
    end
    
    filename = '';
    % Save results if toggle_save is enabled
    if toggle_save
        filename = saveResults(u, ul, y, yl, descriptions, systype, mode, T_sim);
    end

    % Final output formatting
    u = cell2mat(u(:));
    ul = cell2mat(ul(:));
end

%% Helper Functions

function validateInputs(mode, vals)
    % Validate user inputs for mode and vals.
   
    valid_modes = {'r', 'nvstini', 'nt', 'constraints', 'constr', 'mixed'};    
    if ~ischar(mode) || ~ismember(lower(mode), valid_modes)
        error("Invalid mode. Supported modes are: %s", strjoin(valid_modes, ', '));
    end

    required_fields = {'r', 'NvsTini', 'constraints', 'mixed'};
    missing_fields = setdiff(required_fields, fieldnames(vals));
    if ~isempty(missing_fields)
        error("The 'vals' structure is missing required fields: %s", strjoin(missing_fields, ', '));
    end
end

function [values, nruns, param_settings] = resolveMode(mode, vals, max_tries)
    % Resolve mode-specific configurations and parameter settings.
    switch lower(mode)
        case 'r'
            values = vals.r;
            nruns = numel(values);
            param_settings = struct('constraint_scaler', 2, 'T_ini', -1, 'N', -1, 'max_tries', max_tries);
        case {'nvstini', 'nt'}
            values = vals.NvsTini;
            nruns = size(values, 1);
            param_settings = struct('constraint_scaler', 1, 'R', -1, 'max_tries', max_tries);
        case {'constraints', 'constr'}
            values = vals.constraints;
            nruns = numel(values);
            param_settings = struct('T_ini', -1, 'N', -1, 'R', -1, 'max_tries', max_tries);
        case 'mixed'
            values = vals.mixed;
            nruns = size(values.nt, 1) * numel(values.constr) * numel(values.R);
            param_settings = struct('max_tries', max_tries);
        otherwise
            error("Unsupported mode.");
    end
end

function [param_desc, param_values] = buildDescription(mode, param_settings, values, index, systype, T_sim)
    % Build a description and extract parameter values for the given mode.
    switch lower(mode)
        case 'r'
            param_desc = sprintf('ddsfTuner-%s-N%d-Tini-%d-R%d-%s-T%d', mode, ...
                param_settings.N, param_settings.T_ini, values(index), systype, T_sim);
            param_values = struct('N', param_settings.N, 'T_ini', param_settings.T_ini, ...
                'constraint_scaler', param_settings.constraint_scaler, 'R', values(index));
        case {'nvstini', 'nt'}
            param_desc = sprintf('ddsfTuner-%s-N%d-Tini-%d-%s-T%d', mode, ...
                values(index, 2), values(index, 1), systype, T_sim);
            param_values = struct('N', values(index, 2), 'T_ini', values(index, 1), ...
                'constraint_scaler', param_settings.constraint_scaler, 'R', param_settings.R);
        case {'constraints', 'constr'}
            param_desc = sprintf('ddsfTuner-%s-scalingfactor%d-%s-T%d', mode, ...
                values(index), systype, T_sim);
            param_values = struct('N', param_settings.N, 'T_ini', param_settings.T_ini, ...
                'constraint_scaler', values(index), 'R', param_settings.R);
        case 'mixed'
            [i, j, l] = ind2sub([size(values.nt, 1), numel(values.constr), numel(values.R)], index);
            % fprintf('Run %d: i = %d, j = %d, l = %d\n', index, i, j, l);
            param_desc = sprintf('ddsfTuner-%s-Tini%d-N%d-constr_scale%d-R%d-systype-%s-T%d', ...
                mode, values.nt(i, 1), values.nt(i, 2), values.constr(j), values.R(l), systype, T_sim);
            param_values = struct('N', values.nt(i, 2), 'T_ini', values.nt(i, 1), ...
                'constraint_scaler', values.constr(j), 'R', values.R(l));
    end
end

function [u, ul, y, yl, desc] = executeDDSF(systype, T_sim, param_values, desc)
    % Execute the DDSF process and extract logs.
    [~, ~, logs] = runDDSF(systype, T_sim, param_values.N, param_values.T_ini, param_values.constraint_scaler, param_values.R, false);
    u = logs.u;
    ul = logs.ul_t;
    y = logs.y;
    yl = logs.yl;
end

function filename = saveResults(u, ul, y, yl, descriptions, systype, mode, T_sim)
    % Save results to CSV files.
    
    % Define the base output directory
    currentFilePath = mfilename('fullpath');
    currentFolder = fileparts(currentFilePath);
    parentFolder = fileparts(currentFolder);
    base_output_dir = fullfile(parentFolder, 'outputs', 'data'); % Adjusted path construction

    % Ensure the directory exists
    if ~exist(base_output_dir, 'dir')
        mkdir(base_output_dir); % Create the directory if it doesn't exist
    end

    % Construct prefixes for filenames
    prefix_u = fullfile(base_output_dir, sprintf('U-ddsfTuner-systype-%s-mode-%s-T%d', systype, mode, T_sim));
    prefix_y = fullfile(base_output_dir, sprintf('Y-ddsfTuner-systype-%s-mode-%s-T%d', systype, mode, T_sim));

    % Use csvFlexSave to save the data
    uname = csvFlexSave(prefix_u, u, ul, descriptions);
    yname = csvFlexSave(prefix_y, y, yl, descriptions);

    % Return the filenames as a struct
    filename = struct( ...
        'u', uname, ...
        'y', yname ...
    );
end

