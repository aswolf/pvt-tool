% estParamCov - estimate param. covariance from numerical hess of -LogP fun
%
function [pfitcov] = estParamCov(nLogPFun,pfit,opt)
    % To be implemented later
    %opt.pscl
    %opt.credLvl
    %opt.fixFlag
    %optDefault = getEstParamDefaultOpt();

    % For now pscl is unused (eventually it will be taken as an optional
    % argument
    % Use this method for now, since it often works
    [pfitcov,scale,H,scaleFun,pfitUpdate] = calcHessScale(nLogPFun,pfit);
end
