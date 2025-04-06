function r = custom_rank(H)
%   CUSTOM_RANK - Computes the rank of a given matrix H with additional
%                 error handling.

    tolerance = 1e-10; % Small threshold for singular values

    try 
        r = rank(H);
    catch ME
        % Error handling for to sparse matrices
        if contains(ME.message, "SVD does not support sparse matrices")
            % Use svds as a fallback for sparse matrices
            singular_values = svds(H, min(size(H)));
            r = sum(singular_values > tolerance);
        else
            % Re-throw any other errors
            rethrow(ME);
        end
    end
end

