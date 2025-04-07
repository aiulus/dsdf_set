function test_optDDSF()
    % === Load a simple system (with CORA sets)
    systype = 'damper';
    sys = systemsDDSF(systype);
    dims = sys.dims;

    T_ini = sys.config.T_ini;
    N = sys.config.N;
    L = T_ini + N + T_ini;

    %% Setup a dummy lookup structure
    lookup = struct();
    lookup.sys = sys;
    lookup.dims = dims;
    lookup.config = sys.config;
    lookup.IO_params.verbose = false;
    lookup.opt_params = struct( ...
        'R', eye(dims.m), ...
        'regularize', true, ...
        'target_penalty', false, ...
        'constr_type', 'f', ...
        'solver_type', 'o', ...
        'init', true ...
    );

    % Required by gendataDDSF -> inputSignalGenerator
    lookup.data_options = struct( ...
        'datagen_mode', 'uniform', ...
        'scale', 2, ...
        'safe', false ...
    );

    lookup.T_d = 0;

    %% Generate PE input/output data
    [u_d, y_d, ~, ~, ~] = gendataDDSF(lookup);

    %% Construct Hankel matrices
    [H_u, H_y] = hankelDDSF(u_d, y_d, lookup);
    lookup.H = [H_u; H_y];
    lookup.H_u = H_u;
    lookup.H_y = H_y;
    lookup.dims.hankel_cols = size(H_u, 2);

    %% Initial trajectory
    u_ini = u_d(:, 1:T_ini);
    y_ini = y_d(:, 1:T_ini);
    traj_ini = [u_ini; y_ini];

    %% Candidate backup control trajectory
    u_l = u_d(:, 1:N);  % Can be any feasible signal

    %% Run the optimization
    fprintf("🧪 Running optDDSF test...\n");
    [u_opt, y_opt] = optDDSF(lookup, u_l, traj_ini);

    %% === Assertions ===
    assert(~isempty(u_opt), '❌ u_opt is empty');
    assert(~isempty(y_opt), '❌ y_opt is empty');
    assert(all(size(u_opt) == [dims.m, L]), '❌ u_opt size mismatch');
    assert(all(size(y_opt) == [dims.p, L]), '❌ y_opt size mismatch');

    %% Optional: constraint satisfaction (relaxed)
    epsilon = 1e-3;
    U_set = sys.sets.constraints.U;
    Y_set = sys.sets.constraints.Y;

    for k = (T_ini+1):(L-T_ini)
        assert(contains(U_set, u_opt(:, k), epsilon), ...
            sprintf('❌ u_opt(:, %d) violates U constraint', k));
        assert(contains(Y_set, y_opt(:, k), epsilon), ...
            sprintf('❌ y_opt(:, %d) violates Y constraint', k));
    end

    fprintf("✅ optDDSF test passed for '%s'\n", systype);
end
