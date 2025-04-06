function [ranked_scores, ranked_configs] = rankScores(scores, configs)
    if length(scores) ~= length(configs)
        error('Scores and configurations must have the same length.');
    end

    scores_numeric = cell2mat(scores);

    [ranked_scores, sort_indices] = sort(scores_numeric, 'descend');
    
    ranked_configs = configs(sort_indices);
end