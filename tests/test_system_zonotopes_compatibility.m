function test_system_zonotopes_compatibility()
% TEST_SYSTEMSDDSF_ZONOTOPE_COMPATIBILITY
%   Verifies that each system in the DDSF library instantiates:
%     - A valid ss_object (Control System Toolbox)
%     - Valid CORA zonotope fields: X0, U_zono, Y_zono, S_f_zono

    % List of available systems to test
    system_list = {
        'quadrotor', 'damper', 'inverted_pendulum', 'dc_motor', ...
        'cruise_control', 'acc', 'ballNbeam', 'double_pendulum'
    };

    fprintf("=== Running DDSF Zonotope Compatibility Tests ===\n");

    for k = 1:length(system_list)
        sys_name = system_list{k};
        fprintf("\n[✓] Testing system: %s\n", sys_name);

        try
            sys = systemsDDSF(sys_name);

            % --- Check Control System Toolbox Object ---
            assert(isa(sys.ss_object, 'ss'), 'Missing or invalid ss_object');

            % --- Check CORA zonotopes ---
            assert(isa(sys.X0, 'zonotope'), 'Missing or invalid X0 zonotope');
            assert(isa(sys.constraints.U_zono, 'zonotope'), 'Missing or invalid U_zono');
            assert(isa(sys.constraints.Y_zono, 'zonotope'), 'Missing or invalid Y_zono');

            % Optional check: Terminal safe set zonotope
            if isfield(sys, 'S_f_zono')
                assert(isa(sys.S_f_zono, 'zonotope'), 'Invalid S_f_zono zonotope');
            end

            fprintf("    → [PASS] %s\n", sys_name);

        catch ME
            fprintf("    → [FAIL] %s: %s\n", sys_name, ME.message);
        end
    end

    fprintf("\n=== Compatibility Test Complete ===\n");
end
