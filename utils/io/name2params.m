function [descriptions, mode, sysname, sys, dims, T_sim, nruns] = name2params(filename)
    [~, descriptions] = csvFlexRead(filename);
    
    mode = filename2param(filename, 'alg');
    sysname = filename2param(filename, 'sys');
    sys = initializeSys(mode, sysname);
    dims = sys.dims; 
    
    T_sim = filename2param(filename, 't');
    nruns = numel(unique(descriptions));      
end

