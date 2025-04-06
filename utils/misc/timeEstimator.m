function timeEstimator(elapsed, t, T)
    % t -   Current iteration
    % T -   Total #iterations
    avg_per_run = elapsed / t - 1;
    rem = avg_per_run * (T - t +1);
    mins = floor(rem / 60);
    secs = mod(rem, 60);
    if mins >= 60
        hrs = floor(mins /60);
        mins = mod(mins /60);
        fprintf('------------------- Estimated time remaining: %d h, %d m, %.0f s remaining\n', hrs, mins, secs);
        if hrs >= 24
            days = floor(hrs / 24);
            hrs = mod(hrs / 24);
            fprintf('------------------- Estimated time remaining: %d d, %d h, %d m, %.0f s remaining\n', days, hrs, mins, secs);
        end
    end
    fprintf('------------------- Estimated time remaining: %d m, %.0f s remaining\n', mins, secs);
end