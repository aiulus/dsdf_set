function isPE = PEness_check(H)
%   PEness_check Verifies if the Hankel matrix H has full rank
    r = rank(H);
    % r = page_rank(H);
    d = min(size(H));
    if r == d
        isPE = true;
        disp('Hankel matrix has full rank.')
    else
        isPE = false;
        disp("Hankel matrix doesn't have full rank!");
    end
end

