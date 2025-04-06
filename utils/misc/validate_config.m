function validate_config(config, A, C)
    if config.N_p <= config.T_ini
        error("Prediction Horizon (current value: N_p = %d) must be " + ...
            "greater than the length of the initial trajcetory " + ...
            "(current value: T_{ini} = %d)!", sys.config.N_p, sys.config.T_ini);
    end
    lat = system_latency(A, C);
    if lat > config.T_ini
        error("T_ini !>= latency(A,C), but T_ini = %d " + ...
            "and latency(A,C) = %d", sys.config.T_ini, lat);
    end
end