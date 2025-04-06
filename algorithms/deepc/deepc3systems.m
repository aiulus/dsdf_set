function sys = deepc3systems(sys_type)
    discretize = false;
    switch sys_type
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
                'u_min', -inf, ... % Minimum force
                'u_max', inf, ... % Maximum force
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

            A = 1 - (params.damping * params.dt) / params.mass;
            B = params.dt / params.mass;
            C = 1;
            D = 0;

            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
                );

            config = struct( ...
                'T', 41, ... % Window length (default: 41) - This reassigned in the main entry point for DeePC !!
                'T_ini', 5, ... % Initial trajectory length
                'N', 15, ... % Prediction horizon (default: 15)
                's', 2 ... % Sliding length
            );

            opt_params = struct( ...
                        'Q', 100 * eye(size(sys.C, 1)), ... % Output cost matrix (150000)
                        'R', 1e-4 * eye(size(sys.B, 2)) ... % Input cost matrix (0.1)
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
            C = [1 0 0 0;
                 0 0 1 0];
            D = [0;
                 0];

            % Discretize the continuous-time system
            [Ad, Bd, Cd, Dd] = discretize_system(A, B, C, D, params.dt);

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
            C = [1 0];
            D = 0;

            % Discretize the continuous-time system
            [Ad, Bd, Cd, Dd] = discretize_system(A, B, C, D, params.dt);

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
               'dt', 0.1, ... % Sampling time
                'u_min', -100, ... 
                'u_max', 100, ...
                'y_min', -inf, ...
                'y_max', inf, ...
                'x_ini', [0.5;0.5], ... % y_ini = x_ini(1)
                'target', 5,...
                'mass', 1, ...
                'spring_constant', 1, ...
                'damping_coeff', 0.2 ...
                );

            dt = params.dt;
            m = params.mass;
            b = params.damping_coeff;
            k = params.spring_constant;

            % State-space matrices
            A = [1 dt; -k/m*dt 1 - b/m*dt];
            B = [0; dt/m];
            C = [1 0];
            D = 0;    

            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
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
            dt = params.dt;

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

        %% Case 7: Continuous Stirred-Tank Reactor
        case 'cstr'
             params = struct( ...
                'V', 1.0, ... % Reactor volume [m^3]
                'F', 1.0, ... % Volumetric flow rate [m^3/s]
                'k0', 1.0e6, ... % Pre-exponential factor [1/s]
                'E', 50000, ... % Activation energy [J/mol]
                'R', 8.314, ... % Universal gas constant [J/(mol·K)]
                'dH', -50000, ... % Heat of reaction [J/mol]
                'rho', 1000, ... % Density [kg/m^3]
                'Cp', 4.18, ... % Heat capacity [J/(kg·K)]
                'UA', 5000, ... % Heat transfer coefficient [W/K]
                'CA0', 1.0, ... % Inlet concentration [mol/m^3]
                'T0', 300, ... % Inlet temperature [K]
                'Tc', 290, ... % Coolant temperature [K]
                'x_ini', [0.5; 155], ... % Initial state [CA; T]
                'target', [0.5; 310], ... % Desired output [CA; T]                
                'y_min', [-inf;-inf], ... 
                'y_max', [inf;inf], ... 
                'u_min', [-inf;-inf], ... 
                'u_max', [inf;inf] ... 
            );

            A = [   -5  -0.3427; 
                 47.68    2.785];
            B = [    0   1
                   0.3   0];
            C = [0 1
                 1 0];
            D = zeros(2,2);

            % Define the system structure
            sys = struct( ...
                'A', A, ...
                'B', B, ...
                'C', C, ...
                'D', D ...
            );

            % Optimization parameters
            opt_params = struct( ...
                'Q', 1 * eye(size(sys.C, 1)), ... % Output cost matrix
                'R', 0.01 * eye(size(sys.B, 2)) ... % Input cost matrix
            );

            % Run configuration
            config = struct( ...
                'T', 50, ...    % Window length
                'T_ini', 15, ... % Initial trajectory length
                'N', 25, ...  % Prediction horizon
                's', 2 ...      % Sliding length
            );
    end
  
    sys = deepc3_populate_system(sys, params, opt_params, config);

    if discretize
        [sys.A, sys.B, sys.C, sys.D] = discretize_system(A, B, C, D, dt);
    end
end

function sys = deepc3_populate_system(sys, params, opt_params, config)
    sys = constraint_handler(sys, params);
    % Assign config. & optimization parameters
    sys.config = config;
    sys.opt_params = opt_params;
end

