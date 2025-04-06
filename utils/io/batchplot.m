function batchplot(filename)
    [~, descriptions] = csvFlexRead(filename);

    mode = filename2param(filename, 'alg');
    sysname = filename2param(filename, 'sys');
    sys = initializeSys(mode, sysname);
    dims = sys.dims; 

    T_sim = filename2param(filename, 't');
    nruns = numel(unique(descriptions));

    plotData = csv2plotdata(filename, mode, T_sim, dims);

    for i=1:nruns
        configname = generateConfigName(descriptions, i);
        sorted = sortData(plotData, mode, i);
        gridPlot(mode, configname, sys, sorted);
    end
end





















