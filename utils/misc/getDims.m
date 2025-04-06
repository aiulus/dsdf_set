function dims = getDims(A, B, C, D)
    % Assign dimensions
    dims = struct( ...
        'n', size(A, 1), ... % System state dim.
        'm', size(B, 2), ... % Input dim.
        'p', size(C, 1), ... % Output shape
        'q', size(B, 2) + size(C, 1) ... % #I/O variables
    );
end