function sigma_x = normalized_variance(x)
    window_size = 5; % Define window size
    normalized_data = zeros(size(x));

    for t = window_size:length(x)
        window = x(t-window_size+1:t);
        mu_t = mean(window);
        sigma_t = std(window);
        normalized_data(t) = (x(t) - mu_t) / sigma_t;
    end

    sigma_x = var(normalized_data);
end

