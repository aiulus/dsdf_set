function result = vectorFilter(testarr, x_e)
    % Initialize the result as a cell array of the same size as testarr
    result = cell(size(testarr));
    
    % Loop through each entry in x_e
    for i = 1:length(x_e)
        % Check if the entry in x_e is symbolic
        if isa(x_e(i), 'sym')
            % Convert the symbolic expression to a string
            entry_value = char(x_e(i));
            
            % Check for patterns in the symbolic expression
            if contains(entry_value, 'z')
                % If the entry contains 'z', keep the corresponding entry from testarr
                result{i} = testarr(i);
            elseif contains(entry_value, 'k')
                % If the entry contains 'k', substitute k with 1 and evaluate
                syms k
                result{i} = double(subs(x_e(i), k, 1));
            else
                % For other symbolic entries, attempt to evaluate them to a double
                result{i} = double(x_e(i));
            end
        else
            % If the entry is not symbolic, directly convert it to double
            result{i} = double(x_e(i));
        end
    end
    
    % Convert the cell array to a numeric array
    result = cell2mat(result);
end
