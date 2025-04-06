% SYSTEMSDDSF   Initializes predefined system configurations for DDSF.
%
%    SYS = SYSTEMSDDSF(SYS_TYPE, DISCRETIZE) returns the state-space
%    representation and configuration parameters for various predefined
%    dynamic systems. The function provides a convenient way to set up
%    systems for control design, simulation, or data-driven safety filter
%    implementation.
%
%    Inputs:
%       SYS_TYPE   - A string specifying the type of system to initialize.
%                    Available options are:
%                      - 'test': Minimal example system for testing.
%                      - 'example0': Simple state-space model.
%                      - 'cruise_control': Vehicle velocity control model.
%                      - 'quadrotor': 6-DOF quadrotor dynamics.
%                      - 'inverted_pendulum': Cart-pendulum stabilization.
%                      - 'dc_motor': DC motor electrical-mechanical dynamics.
%                      - 'damper': Mass-spring-damper system.
%                      - 'thermostat': Room temperature control.
%                      - 'cstr': Continuous stirred-tank reactor (CSTR).
%       DISCRETIZE - A boolean indicating whether the system should be
%                    discretized using the sampling time (dt) provided in
%                    system parameters. Set to true for discrete-time models.
%
%    Outputs:
%       SYS        - A structure containing:
%                      - A: State-transition matrix (state-space model).
%                      - B: Control input matrix.
%                      - C: Output matrix.
%                      - D: Feedthrough matrix.
%                      - params: System-specific parameters (e.g., mass,
%                        damping, inertia, constraints).
%                      - config: Runtime configuration (e.g., horizon lengths,
%                        window lengths, conservatism).
%
%    See also: DISCRETIZE_SYSTEM, POPULATE_SYSTEM_STRUCT, CONSTRAINT_HANDLER.

function sys = systemsDDSF(sys_type)
    switch sys_type
        %% Example 1: The Quadrotor
        case 'quadrotor'
            discretize = false;
            % System-specific parameters
            params = struct( ...
                'mass', 0.2, ... % Quadrotor mass [kg]
                'g', 9.81, ... % Gravity constant
                'dt', 0.1, ... % Time step for discretization
                'u_min', (1)*(-1)*[1; 0.1; 0.1; 0.1], ... % Minimum force
                'u_max', (1)*[1; 0.1; 0.1; 0.1], ... % Maximum force
                'y_min', (1)*(-1)*[0.2; 0.2; 0.2; 1; 1; 1], ... % Output constraints
                'y_max', (1)*[0.2; 0.2; 0.2; 1; 1; 1], ...  % Output constraints
                'I', repmat(10^(-3), 3, 1), ... % Moment of inertia in x, y, z
                'p', 6, ... % Output dimension (y € R^p)
                'm', 4, ... % Input dimension (u € R^m)
                'n', 12, ... % State dimension (x € R^n)
                'x_ini', zeros(12, 1), ...
                'target', ones(6, 1) ... % TODO: Current value is just a placeholder
                );
    
            run_config = struct( ...
                'T', 214, ... % Data length
                'T_ini', 2, ... % Initial trajectory length
                'N', 20, ... % Prediction horizon
                's', 2 ... % Conservatism
                );
    
            %% State-space Matrices
            % Define state-space matrices as sparse for efficiency
            A_i = [1, 2, 3, 10, 11, 12, 8, 7];
            A_j = [4, 5, 6, 7, 8, 9, 1, 2];
            A_val = [ones(6, 1); params.g; -params.g];
            A = sparse(A_i, A_j, A_val, params.n, params.n);     
    
            B_i = [9, 4, 5, 6];
            B_j = [1, 2, 3, 4];
            B_val = [1/params.mass, 1/params.I(1), 1/params.I(2), 1/params.I(3)];
            B = sparse(B_i, B_j, B_val, params.n, params.m);                              
    
            % Output matrices (position and orientation tracking)
            % Define the indices of x that correspond to y
            indices = [1, 2, 3, 10, 11, 12]; % Indices for ϕ, θ, ψ, x, y, z in x
    
            % Create C as a sparse matrix
            C = sparse(1:length(indices), indices, 1, length(indices), 12);
    
            D = zeros(6, 4);
        
        %% Example 2: Mass Spring damper
        case 'damper'
            discretize = false;
            params = struct( ...
                'dt', 1, ... % Sampling time [s]
                'u_min', 0, ... % U = [0, 100]
                'u_max', 100, ...
                'y_min', -10, ... % Y = [-10, 10]
                'y_max', 10, ...
                'x_ini', [9;2], ... % [vert. displacement, vert. velocity]
                'target', 0,... % vertical displacement [m]
                'mass', 100, ... [kg] 
                'spring_constant', 1, ... [N/m]
                'damping_coeff', 0.2, ... [N*s/m]
                'F', 0 ... [N]
                );
    
            dt = params.dt;
            m = params.mass;
            b = params.damping_coeff;
            k = params.spring_constant;
    
            % State-space matrices
            A = [1,         dt;
                (-(k*dt)/m), (1-(b*dt)/m)];
            B = [0; dt/m];
            C = [1 0]; % Position tracking
            D = 0;
   

            run_config = struct( ...
                'T', 49, ... % Data length
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Conservatism
                );
    
        %% Example 3: Inverted Pendulum
        case 'inverted_pendulum'
            discretize = false;
            params = struct( ...
                'c_mass', 50, ... % Mass of the cart [kg]
                'p_mass', 2, ... % Mass of the pendulum [kg]
                'I', 0.6, ... % Mass moment of inertia of the pendulum [kg.m^2]
                'l', 3, ... % length of the pendulum [m]
                'g', 9.81, ... % Gravity constant [m/s^2]
                'b', 0.1, ... % Friction [N*s/m]
                'dt', 0.1, ... % Time step for discretization
                'y_min', [0;-pi/6], ... % Positional constraint
                'y_max', [1.5;pi/6], ... % Positional constraint
                'u_min', -100, ... % Minimum force
                'u_max', 100, ... % Maximum force
                'target', [1.45, NaN], ... % Desired output
                'x_ini', [0.5; 0; 0.087; 0], ... % Initial state [x, x_dot, theta, theta_dot]
                'state_name', {"Linear Position, Linear Velocity, Angular Position, Angular Velocity"}, ...
                'input_name', {"Force"}); % Initial velocity [m/s]
    
            M = params.c_mass;
            m = params.p_mass;
            I = params.I;
            l = params.l;
            b = params.b;
            g = params.g;
    
            % Compute the state-space matrices
    
            p = I*(M+m)+M*m*l^2; % denominator for the A and B matrices
    
            Ac = [0      1              0           0;
                0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
                0      0              0           1;
                0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
            Bc = [     0;
                (I+m*l^2)/p;
                0;
                m*l/p];
            % Discretize the continuous-time system
            [A, B] = simple_discretize(Ac, Bc, params.dt);
            C = [1 0 0 0;
                0 0 1 0];
            D = [0;
                0];    
           
    
            run_config = struct( ...
                'T', 49, ... % Data length
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Conservatism
                );
    
        %% Example 4: DC Motor
        case 'dc_motor'
            discretize = true;
            params = struct( ...
                'J' , 0.01, ... % Inertia
                'b', 0.1, ... % Damping coefficient
                'K', 0.01, ... % Motor constant
                'R', 1, ... % Resistance
                'L', 0.5, ... % Inductance
                'dt', 0.1, ... % Sampling time
                'u_min', 0, ... % Voltage limits
                'u_max', 24, ... % Voltage limits
                'y_min', 0, ... % Speed limits
                'y_max', 300, ... % Speed limits
                'x_ini', [1; 1], ... % y_ini = x_ini(1)
                'target', 10 ...
                );
    
            b = params.b;
            J = params.J;
            K = params.K;
            R = params.R;
            L = params.L;
    
            A = [-b/J K/J; -K/L -R/L];
            B = [0; 1/L];
            C = [1 0];
            D = 0;
    
            run_config = struct( ...
                'T', 49, ... % Data length
                'T_ini', 2, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Conservatism
                );
    
            %% Example 5: Cruise Control
        case 'cruise_control'
            discretize = false;
            % System-specific parameters
            params = struct( ...
                'mass', 1000, ... % Vehicle mass [kg]
                'damping', 50, ... % Damping coefficient [N*s/m]
                'dt', 0.1, ... % Sampling rate for discetization [s]
                'u_min', 0, ... % Minimum force
                'u_max', 20, ... % Maximum force
                'y_min', -200, ... % Output constraint
                'y_max', 200, ... % Output constraint
                'target', 0, ... % Reference velocity [m/s]
                'slack', 1e-2, ... % For relaxation
                'x_ini', 0, ... % Currently not used
                'state_name', {"Velocity"}, ...
                'input_name', {"Force"}); % Initial velocity [m/s]
    
            A = 1 - (params.damping * params.dt) / params.mass;
            B = params.dt / params.mass;
            C = 1;
            D = 0;
    
            run_config = struct( ...
                'T', 150, ... % Data length
                'T_ini', 10, ... % Initial trajectory length
                'N', 30, ... % Prediction horizon
                's', 1 ... % Conservatism; cannot exceed dims.m in the way this is used in the current implementation
                );
    
            %% Example 6: Adaptive Cruise Control with Time-Delay
        case 'acc'
            discretize = true;
            params = struct( ...
                'mc', 1650, ... % Follower car mass [kg]
                'vl', 20, ... % Lead car velocity [m/s]
                'x_ini', 0.1, ... % Initial distance [km]
                'target', 0.2, ... % Target distance [km]
                'u_min', 0, ... % Control input
                'u_max', 2000, ...  % boundaries
                'y_min', 0.2, ... % Distance variation
                'y_max', 100, ...  % boundaries
                'dt', 0.2, ... % Sampling time [s]
                'Td', 3 ... % Time delay / [dt]
                );
    
            run_config = struct( ...
                'T', 49, ... % Data length
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Conservatism
                );
    
            A = [0 1; 0 0];
            B = [0; (1/params.mc)];
            C = [1 0];
            D = 0;
    
            %% Example 7: Ball & Beam
        case 'ballNbeam'
            discretize = true;
            params = struct( ...
                'm', 0.11, ... % Mass of the ball [kg]
                'R', 0.015, ... % Radius of the ball [m]
                'd', 0.03, ... % Lever arm offset [m]
                'L', 1, ... % Length of the beam [m]
                'J', 9.99e-6, ... % Ball's moment of inertia [kg*m^2]
                'g', 9.8, ... % Gravitational constant [m/s^2]
                'x_ini', 0, ... % Initial ball position
                'target', 0.5, ... % Desired ball position
                'u_min', -10, ... % Minimum gear angle
                'u_max', 10, ...  % Maximum gear angle
                'y_min', 0, ...
                'y_max', 1 ... % Must be the same as L
                );
    
            b21 = - (params.m * params.g * params.d) /(params.L * ...
                (params.m + (params.J / (params.R^2))));
    
            A = [0 1; 0 0];
            B = [0; b21];
            C = [1 0];
            D = 0;

            run_config = struct( ...
                'T', 490, ... % Data length
                'T_ini', 1, ... % Initial trajectory length
                'N', 5, ... % Prediction horizon
                's', 2 ... % Conservatism
                );
    
        %% Example 8: Double Pendulum
        case 'double_pendulum'
            discretize = true;
            % System-specific parameters
            params = struct( ...
                'M', 100, ... % Base mass [kg]
                'M1', 10, ... % Pendulum 1 mass [kg]
                'M2', 10, ... % Pendulum 2 mass [kg]
                'L1', 2, ... % Length of pendulum 1 [m]
                'L2', 1, ... % Length of pendulum 2 [m]
                'g', 9.81, ... % Gravitational acceleration [m/s^2]
                'dt', 0.1, ... % Sampling time [s]
                'x_ini', [-5; zeros(5, 1)], ... % Initial state (position, velocity)
                'target', zeros(6,1), ...
                'u_min', -inf, ...
                'u_max', inf, ...
                'y_min', -inf, ...
                'y_max', inf ...
                );
    
            % Extract parameters for clarity
            M = params.M;
            M1 = params.M1;
            M2 = params.M2;
            L1 = params.L1;
            L2 = params.L2;
            g = params.g;
    
            % State-space matrices
            A = sparse([
                0, 0, 0, 1, 0, 0;
                0, 0, 0, 0, 1, 0;
                0, 0, 0, 0, 0, 1;
                0, -g*(L1*(M1+M2)+L2*M2)/(L1*M), -g*(L2*M2)/(L1*M), 0, 0, 0;
                0, g*(L1*M1*(M + M1 +M2)+L2*M2*(M+M1))/(M*M1*L1^2), g*M2*(-L1*M+L2*(M+M1))/(M*M1*L1^2), 0, 0, 0;
                0, -g*M2/(L1*M1), g*(L1*(M1+M2)-L2*M2)/(L1*L2*M1), 0, 0, 0
                ]);
    
            B = sparse([
                0;
                0;
                0;
                1/M1;
                -1/(L1*M);
                0
                ]);
    
            C = sparse([
                1, 0, 0, 0, 0, 0;
                1, L1, 0, 0, 0, 0;
                1, L1+L2, L2, 0, 0, 0
                ]);
    
            D = sparse(zeros(3, 1));
    
            % Runtime configuration
            run_config = struct( ...
                'T', 300, ... % Data length
                'T_ini', 5, ... % Initial trajectory length
                'N', 20, ... % Prediction horizon
                's', 2 ... % Conservatism
                );                  
    end

    if discretize
        [A, B, C, D] = discretize_system(A, B, C, D, params.dt);
    end
    
    sys = populate_system_struct(A, B, C, D, params); % Collect all system properties in a single object
    
    sys = constraint_handler(sys, params); % Parse constraints
   
    sys.config = validate_config(run_config, A, C);  % Perform checks on adherence to assumptions

    sys.S_f = setEquilibriaDDSF(sys); % Populate the terminal safe set
end
    
function config = validate_config(config, A, C)
    if config.N <= config.T_ini
        error("Prediction Horizon (current value: N = %d) must be " + ...
            "greater than the length of the initial trajcetory " + ...
            "(current value: T_{ini} = %d)!", sys.config.N, sys.config.T_ini);
    end
    lat = system_latency(A, C);
    if lat > config.T_ini
        error("T_ini !>= latency(A,C), but T_ini = %d " + ...
            "and latency(A,C) = %d", sys.config.T_ini, lat);
    end
    min_length = config.N + 2*config.T_ini;
    config.T = (config.T < min_length)*min_length + (config.T >= min_length)*config.T;
end
    

