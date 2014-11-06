% fitErrModPVT - get ErrModel params, ensuring that errorbars match residuals
% pErrMod = [adjVerr, adjTerr]
function [pfitErrMod, pfitcovErrMod, nLogPFun] = fitErrModPVT(pinitErrMod,...
        priorErrMod,priorcovErrMod,Presid,PmarkDerivs,PsampDerivs,...
        VmarkErr,VsampErr,TErr,measGrpID,opt)

    uniqMeasGrpID = unique(measGrpID);
    NMeasGrp = length(uniqMeasGrpID);

    VVTErr = [VmarkErr VsampErr TErr];


    % Construct derivative matrix for residual Pressure values
    %   NOTE: Temp deriv effects tend to cancel out Temp errors due to
    %   difference between sample and maker pressures.
    PresidDerivs = [PmarkDerivs(:,1) PsampDerivs(:,1) ...
        (PmarkDerivs(:,2)-PsampDerivs(:,2))];
    %Vmark, Vsamp, Terr
    %pErrMod = [adjVerr, adjTerr]
    linTransM = [1 0 0; 1 0 0; 0 1 0];

    [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrMod,linTransM,Presid,PresidDerivs,VVTErr,measGrpID,opt);
    [pfitcovErrMod] = estParamCov(nLogPFun,pfitErrMod,[],opt);
end
