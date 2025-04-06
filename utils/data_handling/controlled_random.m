function u = controlled_random(lookup, lb, ub, m, length)
    if ~isfield(lookup.data_options, 'temporal_corr')
        temporal_corr = 0.8; % Default temporal correlation
    else
        temporal_corr = lookup.data_options.temporal_corr;
    end
    u = zeros(m, length);
    for i = 1:m
        u(i, 1) = lb(i) + (ub(i) - lb(i)) * rand; % Initialize first value
        for t = 2:length
            delta = (1 - temporal_corr) * (ub(i) - lb(i)) * rand; % Random deviation
            u(i, t) = min(ub(i), max(lb(i), temporal_corr * u(i, t-1) + delta));
        end
    end
end