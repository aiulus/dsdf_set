function [u_d, y_d, x_d, u, y] = generate_data(sys, T)
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;

    % Generate a persistently exciting, pseudo-random control input
    PE_input = idinput([T, size(B,2)], 'prbs', [0, 1], [-1,1]).'; % Generates a pseudo-random binary signal
    % PE_input = idinput([T, size(B,2)], 'rgs', [0, 1], [-1,1]).';
    %PE_input = idinput([T, sys.params.m], 'prbs', [0, 1], [-1,1]).'; 

    %disp("PE_input: "); disp(PE_input); % DEBUG STATEMENT

    % Initialize input-output storage
    u_d = zeros(sys.params.m, T);
    y_d = zeros(sys.params.p, T);
    x_d = zeros(sys.params.n, T); % Initial state

    disp("(generate_data) u_d size: "); disp(size(u_d)); % DEBUG STATEMENT
    disp("(generate_data) y_d size: "); disp(size(y_d)); % DEBUG STATEMENT

    % Generate data by simulating the system on random inputs for L steps
    for i = 1:T
        u_d(:, i) = PE_input(:, i);
        x_d(:, i + 1) = A * x_d(:, i) + B * u_d(:, i);
        y_d(:, i) = C * x_d(:, i) + D * u_d(:, i);
    end
    % Flatten the control inputs and outputs
    u = reshape(u_d, [], 1); % Reshapes into (T * sys.params.m) x 1
    y = reshape(y_d, [], 1); % Reshapes into (T * sys.params.p) x 1

    disp("(generate_data) u size: "); disp(size(u)); % DEBUG STATEMENT
    disp("(generate_data) y size: "); disp(size(y)); % DEBUG STATEMENT
end

