function validateInputs(mode, vals)
    % Validates the user inputs for mode and vals.

    % Check mode validity
    valid_modes = {'NvsTini', 'QvsR', 'mixed', 'nt', 'qr'};
    if ~ismember(mode, valid_modes)
        error('Invalid mode "%s". Supported modes are: %s.', mode, strjoin(valid_modes, ', '));
    end

    % Check vals structure
    required_fields = {'NvsTini', 'QvsR', 'mixed'};
    for field = required_fields
        if ~isfield(vals, field{1})
            error('The vals structure is missing the required field "%s".', field{1});
        end
    end
end