function PVTeval = fitErrModPVTeval(PVTeval,errModInitList)


    if(isempty(errModInitList))
        errModInitList = PVTeval.errModList;
    end
    errModPriorList = PVTeval.errModPriorList;
    errModCovPriorList = PVTeval.errModCovPriorList;

    PVTdataList = PVTeval.PVTdataList;
    % NEED to figure out how multiple datasets are handled for error model
    assert(numel(PVTdataList)==1,...
        'Multiple datasets in PVTdataList not yet implemented');
    PVTdata = PVTdataList(1);
    errMode = PVTdata.errMode;

    P = PVTdata.Pmark;
    V = PVTdata.V;
    T = PVTdata.T;
    measGrpID = PVTdata.measGrpID;

    % Update initial total P errobars (given errormodel)
    Psamp = PVTeval.Psamp;
    PsampDerivs = PVTeval.PsampDerivs;
    PVTdata = updatePVTdata(PVTdata,PsampDerivs,errModInitList);
    PerrTot = PVTdata.PErrTot;
    Presid = P-Psamp;

    opt = PVTeval.opt;


    switch errMode
        case 'mark'
            VmarkErr = PVTdata.VmarkErr;
            VsampErr = PVTdata.VErr;
            TErr     = PVTdata.TErr;
            PmarkDerivs = PVTdata.PmarkDerivs;
            [pfitErrModList, pfitcovErrModList, nLogPFunList] = ...
                fitErrModPVT(measGrpID,...
                errModInitList,errModPriorList,errModCovPriorList,...
                Presid,PmarkDerivs,PsampDerivs,VmarkErr,VsampErr,TErr,opt);
        case 'std'
            assert(false,['std error mode fitting of error Model params '...
                'not yet implemented'])
        case 'tot'
            assert(false',['There is no need to fit err mod params '...
                'in std error mode. Do not use this function']);
    end


    % Propagate errors to total errors in pressure
    PVTdata = updatePVTdata(PVTdata,PsampDerivs,pfitErrModList);

    %Update list element
    PVTdataList(1) = PVTdata;

    PVTeval.PVTdataList  = PVTdataList;
    PVTeval.errModList = pfitErrModList;
    PVTeval.errModCovList = pfitcovErrModList;
end
