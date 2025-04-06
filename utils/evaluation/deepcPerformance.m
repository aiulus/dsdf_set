function metric = deepcPerformance(sequence, target, tolerance_factor, weights)
    % DEEPCPERFORMANCE    computes a normalized metric for convergence.
    %
    % Inputs:
    %   sequence        - An array (n x m), where each row is the vector at each step.
    %   target          - A target vector of size 1 x m.
    %   tolerance_factor - A scalar multiplier for tolerance, relative to the norm of the target vector.
    %   weights         - A 2-element vector [w1, w2] for weighting convergence speed and final deviation.
    %
    % Outputs:
    %   metric          - A scalar convergence metric (0 to 1, where 1 is ideal convergence).
    %
    % Example usage:
    %   metric = evaluateConvergenceMetric(sequence, target, 0.01, [0.5, 0.5]);

    % Default weights if not provided
    if nargin < 4
        weights = [0.5, 0.5];
    end

    % Compute the target norm and tolerance
    target_norm = norm(target);
    tolerance = tolerance_factor * target_norm; % Tolerance as a fraction of target norm

    % Initialize metrics
    convergence_step = NaN; % First step where norm diff < tolerance
    final_deviation = NaN;  % Final deviation from the target

    % Iterate through the sequence to find the convergence step
    for step = 1:size(sequence, 1)
        current_norm_diff = norm(sequence(step, :) - target);
        if current_norm_diff < tolerance && isnan(convergence_step)
            convergence_step = step;
        end
    end

    % If convergence is never reached, assign the maximum number of steps
    if isnan(convergence_step)
        convergence_step = size(sequence, 1);
    end

    % Compute final deviation from the target
    final_deviation = norm(sequence(end, :) - target);

    % --- Normalize the metrics ---
    speed_score = 1 - (convergence_step / size(sequence, 1)); % Range [0, 1]

    deviation_score = 1 - (final_deviation / target_norm); % Range [0, 1]

    % Clip deviation_score to [0, 1] to avoid negative values
    deviation_score = max(0, deviation_score);

    % --- Combine the metrics using weights ---
    metric = weights(1) * speed_score + weights(2) * deviation_score;
end
