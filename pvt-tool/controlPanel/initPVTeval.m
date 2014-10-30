function PVTeval = initPVTeval(name,sampEosPrior,opt)
    PVTeval.name         = name;
    PVTeval.fixFlag      = [];
    PVTeval.sampEosPrior = sampEosPrior;
    PVTeval.sampEos      = sampEosPrior;

    PVTeval.Psamp        = [];
    PVTeval.Psampderivs  = [];
    PVTeval.Presid       = [];

    PVTeval.PVTdataList  = [];

    PVTeval.errMod       = [];
    PVTeval.errModCov    = [];

    PVTeval.opt          = opt;
end
