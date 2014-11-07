% initPVTdata - create PVTdata struc
function PVTdata = initPVTdata(name,opt)
    assert(isstr(name),'name must be a string name for this PVT dataset.')

    optDefault = getPVTdataDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    PVTdata.name     = name;
    PVTdata.opt      = opt;

    PVTdata.markLbl   = {};
    PVTdata.markEos   = [];

    PVTdata.measGrpID= [];   

    PVTdata.Pmark    = [];
    PVTdata.Vmark    = [];
    PVTdata.V        = [];
    PVTdata.T        = [];

    PVTdata.PmarkErr = [];
    PVTdata.VmarkErr = [];
    PVTdata.VErr     = [];
    PVTdata.TErr     = [];

    PVTdata.errMode  = '';
    PVTdata.PErrTot  = [];

    PVTdata.PmarkDerivs = [];    
end
