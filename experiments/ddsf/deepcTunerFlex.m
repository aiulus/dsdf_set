function [u, y, descriptions, filename] = deepcTunerFlex(mode, vals, systype, T_sim, toggle_save)
    % Main entry point for the DeepC Tuner.
    %
    % INPUTS:
    %   mode        - String specifying the tuning mode ('NvsTini', 'QvsR', 'mixed').
    %   vals        - Struct containing pre-defined values for parameters.
    %   systype     - String specifying the system type.
    %   T_sim       - Simulation time.
    %   toggle_save - (Optional) Boolean to toggle saving results.
    %
    % OUTPUTS:
    %   u           - Cell array of control inputs for each run.
    %   y           - Cell array of system outputs for each run.
    %   descriptions - Cell array of descriptions for parameter configurations.

    % Default value for toggle_save
    if nargin < 5, toggle_save = true; end    

    % Validate inputs
    validateInputs(mode, vals);

    % Determine configurations based on mode
    [values, nruns, param_struct] = getModeConfig(mode, vals);

    % Initialize outputs
    u = cell(1, nruns);
    y = cell(1, nruns);
    descriptions = cell(1, nruns);

    % Output directory setup
    output_dir = fullfile('..', 'outputs', 'plots', 'deepc_tuner');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Main loop for parameter tuning
    tic;
    for i = 1:nruns
        fprintf('\n ------------------- Trying parameter conf. %d / %d -------------------\n', i, nruns);
        if i > 1
            elapsed = toc;
            timeEstimator(elapsed, i, nruns);
        end

        % Extract parameters for this run
        params = extractParameters(mode, i, values, param_struct);

        % Description generation
        d_i = createDescription(mode, systype, T_sim, params);

        % Execute and handle retries
        [u{i}, y{i}] = executeWithRetries(systype, T_sim, params, d_i);

        % Store description
        descriptions{i} = d_i;
    end

    % Save results if toggled
    if toggle_save
        prefix = sprintf('deepcTuner-mode-%s-systype-%s-T%d', mode, systype, T_sim);
        fullpath = getFullPath(prefix);
        filename = csvFlexSave(fullpath, u, y, descriptions);
    end

    % Convert outputs to matrices
    u = cell2mat(u(:));
    y = cell2mat(y(:));
end

function [u_i, y_i] = executeWithRetries(systype, T_sim, params, d_i)
    % Executes the simulation with retries in case of errors.
    success = false;
    max_tries = 5;
    for k = 1:max_tries
        try
            logs = runParamDPC(systype, params.Q, params.R, params.T_ini, params.N, T_sim);
            u_i = logs.u;
            y_i = logs.y;
            success = true;
            return;
        catch ME
            fprintf(['\n Attempt to run runParamDPC.m (conf.: %s) failed at: %s\n ' ...
                     'Message: = %s\n. Trying again...\n'], d_i, ME.stack(1).name, ME.message);
        end
    end
    if ~success
        fprintf('\n Failed after %d attempts for configuration "%s".', max_tries, d_i);
        % Assign fallback values
        u_i = NaN;  
        y_i = NaN;  
    end
end

