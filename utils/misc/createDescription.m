function d_i = createDescription(mode, systype, T_sim, params)
    % Generates a description string for the current configuration.

    if strcmp(mode, 'mixed')
        d_i = sprintf('deepcTuner-%s-N%d-Tini%d-Q%.2f-R%.2f-%s-T%d', ...
                      mode, params.N, params.T_ini, params.Q, params.R, systype, T_sim);
    elseif ismember(mode, {'NvsTini', 'nt'})
        d_i = sprintf('deepcTuner-%s-N%d-Tini%d-%s-T%d', mode, params.N, params.T_ini, systype, T_sim);
    else
        d_i = sprintf('deepcTuner-%s-Q%.2f-R%.2f-%s-T%d', mode, params.Q, params.R, systype, T_sim);
    end
end