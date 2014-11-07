function PVTeval = initPVTeval(name,PVTdataList,sampEosPrior,opt)

    assert(numel(PVTdataList)==1,...
        'multiple data sets in PVTdataList not yet implemented');

    PVTeval.name         = name;
    PVTeval.fixFlag      = [];

    PVTeval.PVTdataList  = PVTdataList;

    PVTeval.sampEosPrior = sampEosPrior;
    PVTeval.sampEosInit  = [];
    PVTeval.sampEosFit   = [];
    PVTeval.nLogPFun     = [];

    PVTeval.Psamp        = [];
    PVTeval.PsampDerivs  = [];

    PVTeval.errMod       = [];
    PVTeval.errModCov    = [];

    PVTeval.opt          = opt;

    [errMod, errModCov] = initErrMod(PVTdataList);
    PVTeval.errMod    = errMod;
    PVTeval.errModCov = errModCov;
end
function [errMod,errModCov] = initErrMod(PVTdataList)
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
    errMod = zeros(NMeasGrp,2);
    errModCredWid = [optPVTdat.errModVmarkPrior optPVTdat.errModVmarkPrior optPVTdat.errModVmarkPrior];
    errModCov = diag(reshape(repmat(errModCredWid.^2,NMeasGrp,1),1,[]));
end
%function [errMod, errModCov] = initErrMod(PVTdataList)
%    Ndat = length(PVTdataList);
%    errMod = zeros(1,2*Ndat);
%
%    errModCov = diag(Inf*ones(size(errMod)));
%end
