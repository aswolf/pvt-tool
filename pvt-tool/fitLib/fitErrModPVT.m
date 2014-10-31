function [pfitErrMod, pfitcovErrMod, nLogPFun] = fitErrModPVT(pinitErrMod,...
        priorErrMod,priorcovErrMod,Presid,PmarkDerivs,PsampDerivs,...
        VmarkErr,VsampErr,TErr,measGrpID,opt)

    uniqMeasGrpID = unique(measGrpID);
    NMeasGrp = length(uniqMeasGrpID);
    assert(length(pinitErrMod) == 3*NMeasGrp, ...
        'errMod must have 3 params for each unique measGrpID value');

    % Construct derivative matrix for residual Pressure values
    %   NOTE: Temp deriv effects tend to cancel out Temp errors due to
    %   difference between sample and maker pressures.
    PresidDerivs = [PmarkDerivs(:,1) PsampDerivs(:,1) ...
        (PmarkDerivs(:,2)-PsampDerivs(:,2))];
    %Vmark, Vsamp, Terr
    VVTErr = [VmarkErr,VsampErr,TErr];
    [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrMod,Presid,PresidDerivs,VVTErr,measGrpID,opt);
    [pfitcovErrMod] = estParamCov(nLogPFun,pfitErrMod,[],opt);
end
