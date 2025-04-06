function H = svdHankel(H, r)
    % SVD
    [W, Sigma, ~] = svd(H, 'econ');
    % Truncate
    W1 = W(:, 1:r);
    Sigma1 = Sigma(1:r, 1:r);
    H = W1 * Sigma1;
end


