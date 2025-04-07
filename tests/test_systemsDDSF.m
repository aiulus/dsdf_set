function test_systemsDDSF()
    %import cora.*

    disp('📦 Running DDSF system tests with CORA integration...');

    % List of supported systems
    sys_list = {'quadrotor', 'damper', 'inverted_pendulum', ...
                'dc_motor', 'cruise_control', 'ballNbeam'};

    for i = 1:length(sys_list)
        sys_type = sys_list{i};
        fprintf('🔍 Testing system: %s\n', sys_type);

        % Call system constructor
        sys = systemsDDSF(sys_type);

        % === Validate Dimensions ===
        n = sys.dims.n;
        m = sys.dims.m;
        p = sys.dims.p;

        assert(isa(sys.sets.X0, 'zonotope'), '❌ X0 not zonotope');
        assert(size(sys.sets.X0.center, 1) == n, '❌ X0 dimension mismatch');

        assert(isa(sys.sets.constraints.U, 'zonotope'), '❌ U not interval');
        assert(length(center(sys.sets.constraints.U)) == m, '❌ U dimension mismatch');

        assert(isa(sys.sets.constraints.Y, 'zonotope'), '❌ Y not interval');
        assert(length(center(sys.sets.constraints.Y)) == p, '❌ Y dimension mismatch');

        % === Terminal Safe Set Test ===
        if isfield(sys.sets, "S_f_zono") && ~representsa(sys.sets.S_f_zono, 'emptySet')
            assert(isa(sys.sets.S_f_zono, 'zonotope'), '❌ S_f_zono not zonotope');
            assert(size(sys.sets.S_f_zono.center, 1) == n, '❌ S_f_zono dimension mismatch');
        else
            warning("⚠️ No terminal zonotope (S_f_zono) computed for %s", sys_type);
        end
    end

    disp('✅ All CORA system integration tests passed.');
end
