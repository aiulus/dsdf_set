function Z = safeZonotope(center, generators, name)
    % safeZonotope - Safe wrapper for CORA.zonotope() creation
    % 
    % Inputs:
    %   center     - vector, zonotope center
    %   generators - matrix, generator matrix
    %   name       - optional string identifier for better warnings
    %
    % Output:
    %   Z - zonotope object if successful, [] otherwise

    if nargin < 3
        name = '';
    else
        name = [' (' name ')'];
    end

    % Validate presence
    if isempty(center) || isempty(generators)
        warning(['Zonotope creation skipped%s: Empty center or generators.'], name);
        Z = [];
        return;
    end

    % Enforce numeric format
    try
        center = double(center(:)); % ensure column vector
        generators = double(generators);
    catch
        warning(['Zonotope creation skipped%s: Cannot convert to numeric.'], name);
        Z = [];
        return;
    end

    % Check for NaN/Inf
    if any(isnan(center(:))) || any(isinf(center(:))) || ...
       any(isnan(generators(:))) || any(isinf(generators(:)))
        warning(['Zonotope creation skipped%s: NaN or Inf in data.'], name);
        Z = [];
        return;
    end

    % Check dimensional consistency
    if size(center, 1) ~= size(generators, 1)
        warning(['Zonotope creation skipped%s: Dim mismatch (center vs. generators).'], name);
        Z = [];
        return;
    end

    % Construct zonotope
    Z = zonotope(center, generators);
end
