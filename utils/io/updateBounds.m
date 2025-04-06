function V_scaled = updateBounds(V, factor)
    V_min = V(:, 1); V_max = V(:, 2);
    factor = (factor - 1) / 2;
    V_scaled_min = V_min - factor .* abs(V_min);
    V_scaled_max = V_max + factor .* abs(V_max);
    V_scaled = [V_scaled_min, V_scaled_max];
end