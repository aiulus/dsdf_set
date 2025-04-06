function u = white_noise(lb, ub, m, length)
    u = zeros(m, length);
    for i = 1:m
        u(i, :) = lb(i) + (ub(i) - lb(i)) .* randn(1, length);
    end
end