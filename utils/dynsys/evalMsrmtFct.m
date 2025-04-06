function y_current = evalMsrmtFct(Gh, x_current, u_current)
    p = numel(Gh);
    y_current = zeros(p, 1);
    for j=1:p
        gj = Gh{j};
        gj_val = gj(x_current, u_current);
        y_current(j) = gj_val;
    end
end