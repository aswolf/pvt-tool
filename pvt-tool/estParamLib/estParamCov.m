% estParamCov - estimate param. covariance from numerical hess of -LogP fun
%
%  Return NaNs as indicator of fixed parameter
function [pfitcov] = estParamCov(nLogPFun,pfit,fixFlag,opt)
    if(isempty(fixFlag))
        fixFlag = zeros(size(pfit));
    end

    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);
    % To be implemented later
    %opt.pscl
    %opt.credLvl
    %opt.fixFlag
    %optDefault = getEstParamDefaultOpt();

    % For now pscl is unused (eventually it will be taken as an optional
    % argument
    % Use this method for now, since it often works
    pFree = pfit(~fixFlag);
    nLogPFreeFun = @(pFree)(nLogPFun(getAllParams(pFree,pfit,fixFlag)));
    [pfitcovFree] = calcHessScale(nLogPFreeFun,pFree);

    pfitcov = zeros(length(pfit));
    [pfit] = getAllParams(pFree,pfit,fixFlag);

    unfitPriorCov = diag(NaN*ones(size(fixFlag)));
    [pfitcov] = updateFreeCov(pfitcovFree,unfitPriorCov,fixFlag);
    indFixDiag = sub2ind(size(pfitcov),find(fixFlag),find(fixFlag));
    pfitcov(indFixDiag) = NaN;
end
