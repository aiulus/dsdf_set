function gridPlot(mode, configname, sys, varargin)
    switch lower(mode)
        case {'u-ddsf', 'ddsf'}
            gridPlotDDSF(lower(mode), configname, sys, varargin{:});
        case 'y-ddsf'
            gridPlotDDSF(lower(mode), configname, sys, varargin{:});
        case 'deepc'
            gridPlotDeePC(configname, sys, varargin{:});
        otherwise
            error('Unsupported mode: %s', mode);
    end
end