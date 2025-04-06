function H = regHankelDDSF(H_u, H_y)
    ru = size(H_u, 1); ry = size(H_y, 1);
    r = min(ru, ry);

    U_svd = svdHankel(H_u, r);
    Y_svd = svdHankel(H_y, r);

    % This won't work for 
    H = [U_svd; Y_svd];
end

