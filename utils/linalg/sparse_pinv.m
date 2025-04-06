function pinv = sparse_pinv(A, k)
    % Computes a pseudo-inverse for sparce matrices
    %
    %   INPUTS:
    %       k - # Singular values to compute (e.g. rank of the problem)

    [U, S, V] = svds(A, k);
    pinv = V * diag(1 ./ diag(S)) * U';
    %pinv = sparse(pinv);
end

