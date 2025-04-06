function sys = populate_system(sys, params, opt_params, config)
    sys = constraint_handler(sys, params);
    sys.config = config;
    sys.opt_params = opt_params;
end