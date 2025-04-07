function test_genDataExtended()
    disp("🔬 Testing genDataExtended...");

    sys = systemsDDSF('damper'); % Simple system
    dims = sys.dims;

    lookup = struct( ...
        'sys', sys, ...
        'dims', dims, ...
        'config', sys.config, ...
        'opt_params', struct('init', true), ...
        'T_d', 0 ...
    );
    lookup.data_options.scale = 1;
    lookup.data_options.datagen_mode = 'uniform';
    lookup.data_options.safe = 1;
    modes = {'single', 'sample', 'set'};

    figure;
    tiledlayout(3, 1);

    for i = 1:numel(modes)
        lookup.uncertainty_options.x0_mode = modes{i};
        data = genDataExtended(lookup);

        nexttile;
        plotDDSFTrajectory(data);
        title(sprintf("DDSF Data - Mode: %s", modes{i}));
    end

    sgtitle("Comparative DDSF Data Generation Modes");
end
