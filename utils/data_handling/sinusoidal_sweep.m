function u = sinusoidal_sweep(lookup, lb, ub, m, length)
    % 4. Sinusoidal Signals with Frequency Sweeping
    if ~isfield(lookup.data_options, 'freq_range')
        freq_range = [0.1, 10]; % Default frequency range (Hz)
    else
        freq_range = lookup.data_options.freq_range;
    end
    t = (0:length-1); u = zeros(m, length);
    for i = 1:m
        freq = freq_range(1) + (freq_range(2) - freq_range(1)) * (t / length);
        u(i, :) = lb(i) + (ub(i) - lb(i)) .* sin(2 * pi * freq .* t);
    end
end