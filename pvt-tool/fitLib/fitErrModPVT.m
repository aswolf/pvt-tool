% fitErrModPVT - get ErrModel params, ensuring that errorbars match residuals
% pErrMod = [adjVerr, adjTerr]
function [pfitErrModList, pfitcovErrModList, nLogPFunList] = fitErrModPVT(measGrpID,...
        pinitErrModList,priorErrModList,priorcovErrModList,...
        Presid,PmarkDerivs,PsampDerivs,VmarkErr,VsampErr,TErr,opt)
    if(isempty(measGrpID))
        measGrpID = '1';
    end
    if(numel(measGrpID) == 1)
        measGrpID = cellstr(num2str(ones(size(Presid))));
    end
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

    [pfitErrModList nLogPFunList] = fitErrModResidMultiGrp(measGrpID,...
        pinitErrModList,priorErrModList,priorcovErrModList,linTransM,...
        Presid,PresidDerivs,VVTErr,opt);
    %[pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
    %    priorErrMod,priorcovErrMod,linTransM,Presid,PresidDerivs,VVTErr,measGrpID,opt);
    %
    if(~iscell(nLogPFunList))
        [pfitcovErrModList] = estParamCov(nLogPFunList,pfitErrModList,[],opt);
    else
        for(j=1:length(nLogPFunList))
            [jpfitcovErrMod] = estParamCov(nLogPFunList{j},pfitErrModList(j,:),[],opt);
            pfitcovErrModList(j,:,:) = jpfitcovErrMod;
        end
    end
end
