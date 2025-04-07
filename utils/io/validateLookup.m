function validateLookup(lookup)
    schema = lookupSchema();
    fields_required = fieldnames(schema);

    for i = 1:numel(fields_required)
        field = fields_required{i};

        if ~isfield(lookup, field)
            error("Missing required field in `lookup`: '%s'", field);
        end

        expected_type = schema.(field);
        actual_value = lookup.(field);

        % Handle special case for chars (which are strings or char vectors)
        if strcmp(expected_type, 'char') && ~(ischar(actual_value) || isstring(actual_value))
            error("Field '%s' must be a char or string.", field);
        elseif strcmp(expected_type, 'struct') && ~isstruct(actual_value)
            error("Field '%s' must be a struct.", field);
        elseif strcmp(expected_type, 'numeric') && ~isnumeric(actual_value)
            error("Field '%s' must be numeric.", field);
        end
    end

    % Warn about unused fields (optional)
    extra = setdiff(fieldnames(lookup), fields_required);
    if ~isempty(extra)
        warning("`lookup` contains extra fields: %s", strjoin(extra, ", "));
    end
end
