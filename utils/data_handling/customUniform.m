%% TODO: DISCARD
function signal = customUniform(factor, lb, ub, d1, d2)
    ub = (ub > 0) .* (factor .* ub) + (ub < 0) .* ((factor^(-1)) .* ub) + (ub == 0) .* (10^factor);
    lb = (lb < 0) .* (factor .* lb) + (lb > 0) .* ((factor^(-1)) .* lb) + (lb == 0) .* (- 10^factor);
    ub = (ub > 0) .* (factor .* ub) + (ub < 0) .* ((factor^(-1)) .* ub);
    lb = (lb < 0) .* (factor .* lb) + (lb > 0) .* ((factor^(-1)) .* lb);
    signal = lb + (ub - lb) .* rand(d1, d2);
end