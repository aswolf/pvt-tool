function [errMod,errModCov] = initErrModPVT(PVTdataList)
    assert(numel(PVTdataList)==1,'Multiple datasets not yet implimented');
    PVTdata = PVTdataList(1);

    measGrpID = PVTdata.measGrpID;
    Ndat = size(PVTdata.V,1);
    optPVTdat = PVTdata.opt;

    uniqMeasGrpID = unique(measGrpID);
    NMeasGrp = length(uniqMeasGrpID);
    if(length(measGrpID)==1)
        val = measGrpID;
        measGrpID = val*ones(Ndat,1);
    end

    % Set prior values for error model
    errMod = zeros(NMeasGrp,3);
    errModCredWid = [optPVTdat.errModVmarkPrior optPVTdat.errModVmarkPrior optPVTdat.errModVmarkPrior];
    errModCov = diag(reshape(repmat(errModCredWid.^2,NMeasGrp,1),1,[]));
end
