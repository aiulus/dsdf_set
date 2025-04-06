function isPE = deepc_PEness_check(u_d, T_ini, N, sys)    
%   deepc_PEness_check - verifies if the data length is sufficient for
%   persistent excitation as described in the DeePC paper.
%   (m + 1) * (T_ini + N + n_B) - 1

    m = size(u_d, 1); % Input dimension
    n_B = size(sys.A, 1); % Minimum repr. dim.
    required = (m + 1) * (T_ini + N + n_B) - 1; % Minimum required data length

    if size(u_d, 2) < required
        fprintf('Data length T=%d is insufficient. Required: %d\n', size(u_d, 2), required);
        isPE = false;
    else
        disp('Data length satisfies persistency of excitation condition.');
        isPE = true;
    end
end

