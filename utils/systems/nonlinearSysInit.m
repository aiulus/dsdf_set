function sys = nonlinearSysInit(sys_type)
    % Initialize system based on type
    switch sys_type
        %% Logistic Map
        case 'logmap'
            params = struct( ...
                'mu', 1, ...             % Nonlinearity parameter
                'x_ini', 2, ...          % Initial state (scalar for logmap)
                'u_min', -5, ...
                'u_max', 5 ...
            );

            syms x u
            % Dynamics
            f1 = params.mu * x * (1 - x);
            g1 = x;

            Fx = {f1}; % Dynamics functions
            Gx = {g1}; % Measurement functions

            % Use symbolic variables
            statevars = x;
            inputvars = u;

            % Validate parameters
            if params.mu <= 0
                error('Parameter mu must be greater than 0 for the logistic map.');
            end

            % System configuration
            config = struct( ...
                'T', 50, ...              % Time horizon
                'T_ini', 10, ...          % Initial trajectory length
                'N', 15 ...               % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);

        %% Van der Pol Oscillator
        case 'van_der_pol'
            params = struct( ...
                'mu', 1, ...             % Nonlinearity parameter
                'x_ini', [2; 0], ...     % Initial state
                'target', 0, ...
                'u_min', -5, ...
                'u_max', 5 ...
            );

            % Dynamics
            syms x1 x2 u
            f1 = x2;
            f2 = params.mu * (1 - x1^2) * x2 - x1 + u;
            g1 = x1;

            Fx = {f1, f2}; % Dynamics functions
            Gx = {g1};     % Measurement functions

            % Use symbolic variables
            statevars = [x1; x2];
            inputvars = u;

            % System configuration
            config = struct( ...
                'T', 100, ...            % Time horizon
                'T_ini', 20, ...         % Initial trajectory length
                'N', 30 ...              % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);

        %% Nonlinear Pendulum
        case 'nonlinear_pendulum'
            params = struct( ...
                'g', 9.81, ...           % Gravity
                'l', 1, ...              % Length of pendulum
                'r', 2, ...
                'x_ini', [pi/4; 0], ...  % Initial angle and angular velocity
                'target', 0, ...
                'u_min', -2, ...
                'u_max', 2 ...
            );

            % Dynamics
            syms x1 x2 u
            f1 = x2;
            f2 = -(params.g / params.l) * sin(x1) - (params.r) * x2;
            g1 = x1;

            Fx = {f1, f2}; % Dynamics functions
            Gx = {g1};     % Measurement functions            

            % Use symbolic variables
            statevars = [x1; x2];
            inputvars = u;

            % System configuration
            config = struct( ...
                'T', 80, ...             % Time horizon
                'T_ini', 10, ...         % Initial trajectory length
                'N', 20 ...              % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);
         %% Inverted Pendulum on a Cart
        case 'inverted_pendulum'
            params = struct( ...
                'g', 9.81, ...         % Gravity
                'l', 1.0, ...          % Length of pendulum
                'm', 0.5, ...          % Mass of pendulum
                'M', 1.0, ...          % Mass of cart
                'x_ini', [0; pi/6; 0; 0], ... % Initial [cart position; pendulum angle; cart velocity; angular velocity]
                'u_min', -10, ...
                'u_max', 10 ...
            );

            % Dynamics
            syms x theta dx dtheta u
            m = params.m; M = params.M; l = params.l; g = params.g;

            % Equations of motion
            f1 = dx;
            f2 = dtheta;
            f3 = (u + m*l*dtheta^2*sin(theta) - m*g*cos(theta)*sin(theta)) / (M + m - m*cos(theta)^2);
            f4 = (g*sin(theta) - cos(theta)*f3) / l;

            g1 = x; % Cart position (output)
            g2 = theta;

            Fx = {f1, f2, f3, f4}; % Dynamics functions
            Gx = {g1, g2}; % Measurement functions

            statevars = [x; theta; dx; dtheta];
            inputvars = u;

            config = struct( ...
                'T', 50, ...          % Time horizon
                'T_ini', 10, ...      % Initial trajectory length
                'N', 15 ...           % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);

        %% Nonlinear Cruise Control with Time Delay
        case 'nonlinear_cruise'
            params = struct( ...
                'tau', 2, ...         % Time delay
                'a', 0.5, ...         % Acceleration gain
                'b', 0.1, ...         % Nonlinear damping
                'x_ini', 20, ...      % Initial velocity
                'u_min', -2, ...
                'u_max', 2 ...
            );

            syms v u u_tau
            % Dynamics with time delay
            f1 = v + params.a*u_tau - params.b*v^2;
            g1 = v; % Output velocity

            Fx = {f1}; % Dynamics functions
            Gx = {g1}; % Measurement functions

            statevars = v; % Current velocity
            inputvars = [u; u_tau]; % Control input and delayed input

            config = struct( ...
                'T', 100, ...         % Time horizon
                'T_ini', 20, ...      % Initial trajectory length
                'N', 30 ...           % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);

        %% Quadcopter with Coupled Nonlinear Dynamics
        case 'quadcopter'
            params = struct( ...
                'g', 9.81, ...        % Gravity
                'm', 0.5, ...         % Mass of quadcopter
                'x_ini', zeros(12,1), ... % Initial state (6 position + 6 velocity states)
                'u_min', -5, ...
                'u_max', 5 ...
            );

            syms x y z phi theta ppsi dx dy dz dphi dtheta dppsi u1 u2 u3 u4
            g = params.g; m = params.m;

            % Dynamics
            f1 = dx; % Velocity in x
            f2 = dy; % Velocity in y
            f3 = dz; % Velocity in z
            f4 = (u1/m) * (cos(phi)*cos(theta)) - g; % z-direction acceleration
            f5 = u2; % Roll dynamics
            f6 = u3; % Pitch dynamics
            f7 = u4; % Yaw dynamics

            Fx = {f1, f2, f3, f4, f5, f6, f7}; % Dynamics functions
            Gx = {x, y, z}; % Position outputs

            statevars = [x; y; z; phi; theta; ppsi; dx; dy; dz; dphi; dtheta; dppsi];
            inputvars = [u1; u2; u3; u4];

            config = struct( ...
                'T', 100, ...         % Time horizon
                'T_ini', 20, ...      % Initial trajectory length
                'N', 25 ...           % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);

        %% Nonlinear Mass-Spring-Damper
        case 'mass_spring_damper'
            params = struct( ...
                'k', 2, ...           % Spring constant
                'c', 1, ...           % Damping coefficient
                'm', 1, ...           % Mass
                'x_ini', [1; 0], ...  % Initial displacement and velocity
                'u_min', -10, ...
                'u_max', 10 ...
            );

            syms x dx u
            k = params.k; c = params.c; m = params.m;

            % Dynamics
            f1 = dx; % Velocity
            f2 = -(k/m)*x - (c/m)*dx^2 + u/m; % Acceleration
            g1 = x; % Displacement output

            Fx = {f1, f2}; % Dynamics functions
            Gx = {g1}; % Measurement functions

            statevars = [x; dx];
            inputvars = u;

            config = struct( ...
                'T', 60, ...          % Time horizon
                'T_ini', 10, ...      % Initial trajectory length
                'N', 15 ...           % Prediction horizon
            );

            opt_params = struct('Q', 1, 'R', 1);
        %% Another Nonlinear Cruise Control with Delay
        case 'acc'
            % System-specific parameters
            params = struct( ...
                'mass', 1000, ...        % Vehicle mass [kg]
                'drag_coeff', 0.32, ...  % Nonlinear drag coefficient
                'dt', 0.1, ...           % Sampling time [s]
                'time_delay', 1, ...   % Variable delay [s]
                'u_min', 0, ...          % Minimum throttle input
                'u_max', 2000, ...       % Maximum throttle input
                'y_min', 0, ...          % Minimum velocity [m/s]
                'y_max', 40, ...         % Maximum velocity [m/s]
                'x_ini', 5, ...          % Initial speed [m/s]
                'target', 20 ...         % Desired speed [m/s]
            );
        
            % Define symbolic variables for state and input
            syms v u u_tau  % Velocity (v), Throttle input (u), Delayed input (u_tau)
        
            % Nonlinear dynamics including delay
            f1 = (u_tau - params.drag_coeff * v^2) / params.mass; % Dynamics
            g1 = v; % Output (velocity)
        
            Fx = {f1}; % Dynamics function
            Gx = {g1}; % Measurement function
        
            statevars = v; % Velocity as the state variable
            inputvars = [u; u_tau]; % Current and delayed inputs
        
            % Linearized system matrices around the desired speed
            v_eq = params.target; % Equilibrium speed
            A = -2 * params.drag_coeff * v_eq / params.mass;
            B = 1 / params.mass;
            C = 1;
            D = 0;
        
            % Create the linearized continuous-time system
            delayss = ss(A, B, C, D);
        
            % Add delay to the linearized system
            delayss.InputDelay = params.time_delay;
            params.delayss = delayss;

            % Define configuration and optimization parameters
            config = struct( ...
                'T', 50, ...          % Time horizon
                'T_ini', 10, ...      % Initial trajectory length
                'N', 15 ...           % Prediction horizon
            );
        
            opt_params = struct('Q', 1, 'R', 1);
        otherwise
            error('Unknown system type specified.');
    end
    
    % Assign the system-specific parameters
    sys.params = params; sys.config = config; sys.opt_params = opt_params;
    [sys.functions.Fh, sys.functions.Gh] = getHandles(Fx, Gx, statevars, inputvars);
    
    try
        sys.S_f = getEquilibriumNonlinear(Fx, statevars, inputvars);
    catch ME
        error('<< getEquilibriumNonlinear >>\n<< Message: >> %s \n', ME.message);
    end

    x_e_default = sys.S_f.trivial_solution.x_e;
    u_e_default = sys.S_f.trivial_solution.u_e;
    [A, B, C, D] = linearize(Fx, Gx, statevars, inputvars, x_e_default, u_e_default);
    sys.A = A; sys.B = B; sys.C = C; sys.D = D;

    [y_e_triv, y_e] = populateYs(sys.S_f, C, D);
    sys.S_f.trivial_solution.y_e = y_e_triv;
    sys.S_f.symbolic_solution.y_e = y_e;   

    sys.dims = getDims(A, B, C, D);
    sys = constraint_handler(sys, params);
    sys.dt = 0.1; % Time step size [s]
    sys.config.s = 1; % Receding horizon-step size
    sys.target = sys.params.target;
end

%% Helper Functions

function [y_e_triv, y_e] = populateYs(S_f, C, D)
    x_e_triv = S_f.trivial_solution.x_e;
    x_e =  S_f.symbolic_solution.x_e;
    u_e_triv = S_f.trivial_solution.u_e;
    u_e =  S_f.symbolic_solution.u_e;
    y_e_triv = C*x_e_triv + D*u_e_triv;
    y_e = C*x_e + D*u_e;    
end

function [Fh, Gh] = getHandles(Fx, Gx, x, u)
    Fh = cell(size(Fx)); Gh = cell(size(Gx));
    for i=1:numel(Fx)
        fi = Fx{i};
        fi_handle =  matlabFunction(fi, 'Vars', {x, u});
        Fh{i} = fi_handle;
    end
    for j=1:numel(Gx)
        gj = Gx{j};
        gj_handle =  matlabFunction(gj, 'Vars',{x, u});
        Gh{j} = gj_handle;
    end
end

