function basis = matrixImage(A)
    [~, pivotCols] = rref(A);
    basis = A(:, pivotCols);
end

