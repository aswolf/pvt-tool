function PVTeval = addDataPVTeval(PVTeval,PVTdataList)
    assert(isempty(PVTeval.PVTdataList),...
        ['Sequential data processing not yet implemented. '...
        'Data must be added all at once'.]);
    assert(numel(PVTdataList)==1,'Multiple datasets not yet implimented');

    PVTeval.PVTdataList = PVTdataList;
    [errMod, errModCov] = initErrModPVT(PVTdataList);
    PVTeval.errMod    = errMod;
    PVTeval.errModCov = errModCov;
end
