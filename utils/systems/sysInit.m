    % SYSINIT   Initializes predefined system configurations.
    %
    %    SYS = SYSINIT(SYS_TYPE) returns the state-space representation and
    %    configuration parameters for various predefined dynamic systems. The
    %    function is designed to simplify system setup for control design,
    %    simulation, and experimentation.
    %
    %    Inputs:
    %       SYS_TYPE   - A string specifying the type of system to initialize.
    %                    Available options include
    %
    %                   'test'              - A basic example system for testing purposes.
    %                   'example0'          - A simple state-space example with predefined matrices.
    %                   'cruise_control'    - A system modeling a vehicle's velocity control.
    %                   'quadrotor'         - A 6-DOF quadrotor model for position and orientation tracking.
    %                   'inverted_pendulum' - An inverted pendulum on a cart for control dynamics.
    %                   'dc_motor'          - A DC motor system modeling electrical and mechanical dynamics.
    %                   'damper'           - A mass-spring-damper system for vibration analysis.
    %                   'thermostat'        - A single-zone temperature control system.
    %                   'cstr'              - A Continuous Stirred-Tank Reactor for chemical process control.
    %
    %    Outputs:
    %       SYS        - A structure containing:
    %                      - A: State-transition matrix (state-space model).
    %                      - B: Control input matrix.
    %                      - C: Output matrix.
    %                      - D: Feedthrough matrix.
    %                      - params: System-specific parameters (e.g., mass,
    %                        damping, constraints, target states).
    %                      - config: Runtime configuration parameters (e.g.,
    %                        horizon lengths, trajectory lengths, conservatism).
    
function sys = sysInit(sys_type)
    switch sys_type
        case 'test'
            params = struct( ...
                'target', [0, 0, 0], ...
                'x_ini', [1; 1; 1] ...
                );
    
            sys = struct( ...
                'A', [1 0 0; ...
                0 1 0;
                0 0 1], ...
                'B', eye(3,2), ...
                'C', eye(3), ...
                'D', zeros(3, 2) ...
                );
    
            config = struct( ...
                'T', 138, ... % Window length
                'T_ini', 5, ... % Initial trajectory length
                'N', 4, ... % Prediction horizon
                's', 2 ... % Sliding length
                );
    
            opt_params = struct( ...
                'Q', eye(size(sys.C, 1)), ... % Output cost matrix
                'R', eye(size(sys.B, 2)) ... % Control cost matrix
                );
    
            %% Case 1: Example0
        case 'example0'
            params = struct( ...
                'target', [0, NaN, NaN], ...
                'x_ini', [1; 1; 1] ...
                );
    
            sys = struct( ...
                'A', [1 -0.01 0; ...
                0.01 1 0;
                0 0 0.8], ...
                'B', eye(3,2), ...
                'C', eye(3), ...
                'D', zeros(3, 2) ...
                );
    
            config = struct( ...
                'T', 138, ... % Window length
                'T_ini', 5, ... % Initial trajectory length
                'N', 4, ... % Prediction horizon
                's', 2 ... % Sliding length
                );
    
            opt_params = struct( ...
                'Q', eye(size(sys.C, 1)), ... % Output cost matrix
                'R', eye(size(sys.B, 2)) ... % Control cost matrix
                );
    
            %% Case 2: Cruise Control
        case 'cruise_control'
            % System-specific parameters
            params = struct( ...
                'mass', 1000, ... % Vehicle mass [kg]
                'damping', 50, ... % Damping coefficient [N*s/m]
                'dt', 0.1, ... % Sampling rate for discetization [s]
                'u_min', 0, ... % Minimum force
                'u_max', 2000, ... % Maximum force
                'y_min', -inf, ... % Output constraint
                'y_max', inf, ... % Output constraint
                'target', 20, ... % Reference velocity [m/s]
                'slack', 1e-2, ... % For relaxation
                'x_ini', 0, ...
                'state_name', {"Velocity"}, ...
                'input_name', {"Force"}); % Initial velocity [m/s]
    
            % 'target', [1, 2, 3, 4] for dims.p = 4
            %  The constraint then just refers to the decision variable at
            %   y(index) being near target(index)
    
            A = 1 - ((params.damping * params.dt) / params.mass);
            B = params.dt / params.mass;
            C = 1;
            D = 0;
    
            %[A, B, C, D] = discretize_system(A, B, C, D, params.dt);
    
            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
                );
    
            config = struct( ...
                'T', 41, ... % Window length - This reassigned in the main entry point for DeePC !!
                'T_ini', 10, ... % Initial trajectory length
                'N', 30, ... % Prediction horizon (default: 15)
                's', 2 ... % Sliding length
                );
    
            opt_params = struct( ...
                'Q', 1e+5 * eye(size(sys.C, 1)), ... % Output cost matrix (150000)
                'R', 1e-2 * eye(size(sys.B, 2)) ... % Input cost matrix (0.1)
                ); % Optimization parameters
        case 'quadrotor'
            % System-specific parameters
            params = struct( ...
                'mass', 0.2, ... % Quadrotor mass [kg]
                'g', 9.81, ... % Gravity constant
                'dt', 0.1, ... % Time step for discretization
                'u_min', (100)*(-1)*[1; 0.1; 0.1; 0.1], ... % Minimum force
                'u_max', (100)*[1; 0.1; 0.1; 0.1], ... % Maximum force
                'y_min', (100)*(-1)*[0.2; 0.2; 0.2; 1; 1; 1], ... % Output constraints
                'y_max', (100)*[0.2; 0.2; 0.2; 1; 1; 1], ...  % Output constraints
                'I', repmat(10^(-3), 3, 1), ... % Moment of inertia in x, y, z
                'p', 6, ... % Output dimension (y € R^p)
                'm', 4, ... % Input dimension (u € R^m)
                'n', 12, ... % State dimension (x € R^n)
                'x_ini', zeros(12, 1), ...
                'target', ones(6, 1) ... % TODO: Current value is just a placeholder
                );
    
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
            indices = [1, 2, 3, 10, 11, 12]; % Indices for ϕ, θ, ψ, x, y, z in x
    
            % Create C as a sparse matrix
            C = sparse(1:length(indices), indices, 1, length(indices), 12);
    
            D = zeros(6, 4);
    
            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
                );
    
            config = struct( ...
                'T', 214, ... % Window length - This reassigned in the main entry point for DeePC !!
                'T_ini', 2, ... % Initial trajectory length
                'N', 20, ... % Prediction horizon (default: 15)
                's', 2 ... % Sliding length
                );
    
            opt_params = struct( ...
                'Q', 1 * eye(size(sys.C, 1)), ... % Output cost matrix (150000)
                'R', 1 * eye(size(sys.B, 2)) ... % Input cost matrix (0.1)
                ); % Optimization parameters
            %% Case 3: Inverted Pendulum
        case 'inverted_pendulum'
            params = struct( ...
                'c_mass', 50, ... % Mass of the cart [kg]
                'p_mass', 2, ... % Mass of the pendulum [kg]
                'I', 0.6, ... % Mass moment of inertia of the pendulum [kg.m^2]
                'l', 3, ... % length of the pendulum [m]
                'g', 9.81, ... % Gravity constant [m/s^2]
                'b', 0.1, ... % Friction [N*s/m]
                'dt', 0.1, ... % Time step for discretization
                'y_min', [0;-inf], ... % Positional constraint
                'y_max', [1.5;inf], ... % Positional constraint
                'u_min', -inf, ... % Minimum force
                'u_max', inf, ... % Maximum force
                'target', [1.45, NaN], ... % Desired output
                'x_ini', [0.5; 0; 0; 0], ... % Initial state [x, x_dot, theta, theta_dot]
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
    
            A = [0      1              0           0;
                0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
                0      0              0           1;
                0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
            B = [     0;
                (I+m*l^2)/p;
                0;
                m*l/p];
            [Ad, Bd] = simple_discretize(A, B, params.dt);
            Cd = [1 0 0 0;
                0 0 1 0];
            Dd = [0;
                0];
    
            % Discretize the continuous-time system
            % [Ad, Bd, Cd, Dd] = discretize_system(A, B, C, D, params.dt);
    
            sys = struct( ...
                'A', Ad, ...
                'B', Bd, ...
                'C', Cd, ...
                'D', Dd ...
                );
    
            opt_params = struct( ...
                'Q', 150000 * eye(size(sys.C, 1)), ... % Output cost matrix
                'R', 0.1 * eye(size(sys.B, 2)) ... % Input cost matrix
                ); % Optimization parameters
    
            config = struct( ...
                'T', 37, ... % Window length
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Sliding length
                );
    
            %% Case 4: DC Motor
        case 'dc_motor'
            % Parameters
            params = struct( ...
                'J' , 0.01, ... % Inertia
                'b', 0.1, ... % Damping coefficient
                'K', 0.01, ... % Motor constant
                'R', 1, ... % Resistance
                'L', 0.5, ... % Inductance
                'dt', 0.1, ... % Sampling time
                'u_min', -inf, ... % Voltage limits
                'u_max', inf, ... % Voltage limits
                'y_min', -inf, ... % Speed limits
                'y_max', inf, ... % Speed limits
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
            [Ad, Bd] = simple_discretize(A, B, params.dt);
            Cd = [1 0];
            Dd = 0;
    
            % Discretize the continuous-time system
            %[Ad, Bd, Cd, Dd] = discretize_system(A, B, C, D, params.dt);
    
            sys = struct( ...
                'A', Ad, ...
                'B', Bd, ...
                'C', Cd, ...
                'D', Dd ...
                );
    
            opt_params = struct( ...
                'Q', 1 * eye(size(sys.C, 1)), ... % Output cost matrix
                'R', 0.1 * eye(size(sys.B, 2)) ... % Input cost matrix
                ); % Optimization parameters
    
            config = struct( ...
                'T', 20, ... % Window length % Not used
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Sliding length
                );
    
        %% Case 5: Mass Spring damper
        case 'damper'
            params = struct( ...
                'dt', 1, ... % Sampling time [s]
                'u_min', 0, ...
                'u_max', 100, ...
                'y_min', -10, ...
                'y_max', 10, ...
                'x_ini', [9;2], ... % [vert. displacement, vert. velocity]
                'target', 0,... % vertical displacement [m]
                'mass', 100, ... [kg] 
                'spring_constant', 1, ... [N/m]
                'damping_coeff', 0.2, ... [N*s/m]
                'F', 10 ... [N]
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
    
            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
                );
    
            opt_params = struct( ...
                'Q', 100 * eye(size(sys.C, 1)), ... % Output cost matrix
                'R', 0.01 * eye(size(sys.B, 2)) ... % Input cost matrix
                ); % Optimization parameters
    
            config = struct( ...
                'T', 20, ... % Window length % Not used
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon
                's', 2 ... % Sliding length
                );
    
            %% Case 6: Temperature Control
        case 'thermostat'
            % System-specific parameters
            params = struct( ...
                'thermal_capacitance', 500, ... % [J/°C]
                'heat_transfer_coeff', 10, ...  % [W/°C]
                'dt', 0.1, ...                  % Sampling time [s]
                'u_min', 0, ...              % Minimum heating power [W]
                'u_max', 15000, ...               % Maximum heating power [W]
                'y_min', 15, ...                % Minimum room temperature [°C]
                'y_max', 27, ...                % Maximum room temperature [°C]
                'target', -5, ...               % Desired room temperature [°C]
                'x_ini', 10 ...                % Initial room temperature [°C]
                );
    
            C = params.thermal_capacitance;
            h = params.heat_transfer_coeff;
    
            % Continuous-time state-space matrices
            A = -h / C;
            B = 1 / C;
            C = 1;
            D = 0;
    
            % Define the system structure
            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
                );
    
            % Optimization parameters
            opt_params = struct( ...
                'Q', 10 * eye(size(sys.C, 1)), ... % Output cost matrix
                'R', 0.01 * eye(size(sys.B, 2)) ... % Input cost matrix
                );
    
            % Run configuration
            config = struct( ...
                'T', 50, ...    % Window length
                'T_ini', 15, ... % Initial trajectory length
                'N', 25, ...  % Prediction horizon
                's', 2 ...      % Sliding length
                );
    
        %% Case 7: Aircraft pitch
        case 'aircraft_pitch'
            % Dynamics equation for d_alpha with delta included
            k_delta = 0.05; % Example coefficient for elevator influence on lift
    
            alpha = deg2rad(5);       % Angle of attack [rad]
            q = 0.05;                 % Pitch rate [rad/s]
            theta = deg2rad(6);       % Pitch angle [rad]
            delta = deg2rad(3);       % Elevator deflection angle [rad]
            rho = 1.225;              % Air density at sea level [kg/m^3]
            S = 16;                   % Wing area [m^2]
            c = 1.5;                  % Mean aerodynamic chord [m]
            m = 750;                  % Aircraft mass [kg]
            U = 50;                   % Equilibrium flight speed [m/s]
            CT = 0.08;                % Coefficient of thrust
            CD = 0.03;                % Coefficient of drag
            CL = 0.5;                 % Coefficient of lift
            CW = CL;                  % Coefficient of weight (assuming lift equals weight)
            CM = -0.02;               % Coefficient of pitch moment
            gamma = deg2rad(3);       % Flight path angle [rad]
            iyy = 0.05;               % Normalized moment of inertia (unitless)
            
            % Compute Dependent Parameters
            mu = (rho * S * c) / (4 * m);         % [1/m]
            Omega = (2 * U) / c;                  % [1/s]
            sigma = 1 / (1 + mu * CL);            % Unitless
            eta = mu * sigma * CM;                % Unitless
            
            % Build Matrices A and B
            A = [
                % d_alpha
                mu * Omega * sigma * (- (CL + CD)), ...
                mu * Omega * sigma * (1 / (mu - CL)), ...
                mu * Omega * sigma * (- CW * sin(gamma));
                
                % d_q
                (mu * Omega / (2 * iyy)) * (CM - eta * (CL + CD)), ...
                (mu * Omega / (2 * iyy)) * (CM + sigma * CM * (1 - mu * CL)), ...
                0;
                
                % d_theta
                0, Omega, 0
            ];
            
            B = [
                % d_alpha
                k_delta;
                
                % d_q
                (mu * Omega / (2 * iyy)) * (eta * CW * sin(gamma));
                
                % d_theta
                0
            ];
            
            C = [0 0 1]; D = 0;

            % Store all parameters in a struct
            params = struct(...
                'alpha', alpha, ...
                'q', q, ...
                'theta', theta, ...
                'delta', delta, ...
                'rho', rho, ...
                'S', S, ...
                'c', c, ...
                'm', m, ...
                'U', U, ...
                'CT', CT, ...
                'CD', CD, ...
                'CL', CL, ...
                'CW', CW, ...
                'CM', CM, ...
                'gamma', gamma, ...
                'iyy', iyy, ...
                'mu', mu, ...
                'Omega', Omega, ...
                'sigma', sigma, ...
                'eta', eta, ...
                'x_ini', [0; 0; 0], ...
                'target', [NaN; NaN; deg2rad(3)], ...
                'u_min', 0, ...
                'u_max', inf, ...
                'y_min', [-inf; -inf; -30], ...
                'y_max', [inf; inf; 30] ...
            );


            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
                );
    
            % Optimization parameters
            opt_params = struct( ...
                'Q', 1 * eye(size(sys.C, 1)), ... % Output cost matrix
                'R', 1 * eye(size(sys.B, 2)) ... % Input cost matrix
                );
    
            % Run configuration
            config = struct( ...
                'T', 250, ...    % Window length
                'T_ini', 3, ... % Initial trajectory length
                'N', 15, ...  % Prediction horizon
                's', 2 ...      % Sliding length
                );  
    
    end
    sys.params = params;
    sys.dims = getDims(sys.A, sys.B, sys.C, sys.D);
    sys = populate_system(sys, params, opt_params, config);
    sys.S_f = setEquilibriaDDSF(sys);
end
    
  

