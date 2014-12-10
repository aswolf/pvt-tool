function PVTeval = initPVTeval(name,PVTdataList,sampEosPrior,opt)

    assert(numel(PVTdataList)==1,...
        'multiple data sets in PVTdataList not yet implemented');

    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    PVTeval.name         = name;
    PVTeval.fixFlag      = [];

    PVTeval.PVTdataList  = PVTdataList;

    PVTeval.sampEosPrior = sampEosPrior;
    PVTeval.sampEosInit  = [];
    PVTeval.sampEosFit   = [];
    PVTeval.nLogPFun     = [];

    PVTeval.Psamp        = [];
    PVTeval.PsampDerivs  = [];

    PVTeval.errModList   = [];
    PVTeval.errModPriorList   = [];
    PVTeval.errModCovPriorList= [];
    PVTeval.errModCovList= [];

    PVTeval.opt          = opt;

    [errModList, errModCovList] = initErrModList(PVTdataList);
    PVTeval.errModList= errModList;
    PVTeval.errModPriorList = errModList;
    PVTeval.errModCovPriorList = errModCovList;
    PVTeval.errModCovList = NaN*errModCovList;
end
function [errModList,errModCovList] = initErrModList(PVTdataList)
    assert(numel(PVTdataList)==1,'Multiple datasets not yet implimented');
    PVTdata = PVTdataList(1);

    measGrpID = PVTdata.measGrpID;
    assert(iscellstr(measGrpID),'measGrpID must be a cell array of str');
    Ndat = size(PVTdata.V,1);
    optPVTdat = PVTdata.opt;

    uniqMeasGrpID = unique(measGrpID);
    NMeasGrp = length(uniqMeasGrpID);
    if(NMeasGrp==1)
        measGrpID = cellstr(repmat(measGrpID,Ndat,1));
    end

    % Set prior values for error model
    errModList = zeros(NMeasGrp,2);
    errModCredWid = [optPVTdat.errModVmarkPrior optPVTdat.errModTPrior];
    errModCovGrp = diag(errModCredWid.^2);
    errModCovList = zeros(NMeasGrp,size(errModCovGrp,1),size(errModCovGrp,2));
    for(i=1:NMeasGrp)
        errModCovList(i,:,:) = errModCovGrp;
    end
end
