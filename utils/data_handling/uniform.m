function signal = uniform(lb, ub, d1, d2)
    signal = lb + (ub - lb) .* rand(d1, d2);
end