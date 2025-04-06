function [ranked_scores, ranked_configs] = evaluateDeePC(filename)
    %, mode, toggle_plot, toggle_save
    [descriptions, mode, ~, sys, dims, T_sim, nruns] = name2params(filename);

    target = sys.params.target;
    
    plotdata = csv2plotdata(filename, mode, T_sim, dims);

    weights = [0.5, 0.5]; % weights(1) for convergence speed,  weights(2) for final result
    
    evaluation_scores = cell(1, nruns); configs = zeros(2, nruns);

    for i=1:nruns
        sorted_i = sortData(plotdata, mode, i);

        y_i = sorted_i.y; y_i = y_i(:);        
            
        m_i = deepcPerformance(y_i, target, 0.8, weights);

        evaluation_scores{i} = m_i;     
        
    end
    [ranked_scores, ranked_configs] = rankScores(evaluation_scores, descriptions);
end

function param = getParam(mode)
    switch mode
        case 'n'
            param = file2param(filename, 'N');
        case 't'
            param = file2param(filename, 'T_ini');
        case 'r'
            param = file2param(filename, 'R');
        case 'q'
            param = file2param(filename, 'Q');
    end
end

