function testDDSFSetBased()
    clc; close all; 

    % === System setup ===
    systype = 'quadrotor';  
    T_sim = 20;
    T_ini = -1; % Use system default
    N = -1;     % Use system default
    scale_constraints = -1; % Keep original bounds
    R = -1; % Use default R
    toggle_plot = 0;

    % === Test modes: 'single', 'sample', 'set' ===
    for mode = ["single", "sample", "set"]
        fprintf("\n==================== Running DDSF with x0_mode = '%s' ====================\n", mode);
        lookup = populateLookup(systype, T_sim, N, T_ini, scale_constraints, R, toggle_plot, mode);
        %% TODO: See what use was originally intended for this field
        lookup.data_options.scale = 1;
        
        % === Generate data ===
        data = genDataExtended(lookup);
        fprintf("Generated data with mode: %s\n", data.mode);

        % === Plot trajectory with zonotope visualization ===
        %figure('Name', ['Output Trajectory - x0\_mode = ', mode]);
        figure('Name', ['Output Trajectory - x0\_mode = ', char(mode)]);
        plotDDSFTrajectory(data);
        title(['DDSF Output Trajectory (', mode, ' mode)']);

        % === (Optional) Run full DDSF receding horizon test ===
        if strcmp(mode, "single")
            fprintf("Running full DDSF closed-loop rollout...\n");
            [~, ~, logs] = runDDSF(systype, T_sim, N, T_ini, scale_constraints, R, 1);
        end
    end
end
