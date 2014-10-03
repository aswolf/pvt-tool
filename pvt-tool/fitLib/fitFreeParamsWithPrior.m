% fitFreeParamsWithPrior - return best-fit params given prior and nLogP Fun
function [pfit nLogPFun] = fitFreeParamsWithPrior(pinit,fixFlag,...
        prior,priorcov,wtResidFun,doRobustFit)

    nLogPFun = makeNLogPFun(wtResidFun,prior,priorcov,doRobustFit,[]);
    
    % Obtain free parameter subset 
    %   - Update fixFlag to include variables with zero variance prior values
    [pinitFree,priorFree,priorcovFree,fixFlag] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag);

    % Test if any free parameters remain
    if(all(fixFlag==1))
        pfit = pinit;
    end

    % Define wrapper function for nLogP to use only free parameters
    nLogPFreeFun = @(pFree)(nLogPFreeFun(getAllParams(pFree,pinit,fixFlag)));
    pfitFree = fitParams(pinitFree,nLogPFreeFun);

    pfit = pinit;
    pfit(~fixFlag) = pfitFree;
end
