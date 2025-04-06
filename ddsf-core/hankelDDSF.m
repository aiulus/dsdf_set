function [H_u, H_y] = hankelDDSF(u_d, y_d, lookup)
    T_ini = lookup.sys.config.T_ini;
    N = lookup.sys.config.N;
    PE_order = N + 2 * T_ini;

    [~, H_u] = construct_hankel(u_d, PE_order);
    [~, H_y] = construct_hankel(y_d, PE_order);

    full_rank = PEness_check(H_u);
    if ~full_rank
        error(['Persistency of excitation check failed. ' ...
            'Please provide richer input data or adjust T_ini and N.']);
    end
end