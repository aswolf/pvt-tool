% fitFreeParamsWithPrior - return best-fit params given prior and nLogP Fun
function [pfit nLogPFun] = fitFreeParamsWithPrior(pinit,fixFlag,...
        prior,priorcov,wtResidFun,opt)

    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    % If prior unspecified, set to Inf in all parameters
    if(isempty(priorcov))
        priorcov = diag(Inf*ones(size(pinit)));
    end

    % Set diagonal value of cov matrix to NaN to indicate fixed parameters
    indFixDiag = sub2ind(size(priorcov),find(fixFlag==1),find(fixFlag==1));
    priorcov(indFixDiag) = NaN;
    nLogPFun = makeNLogPFun(wtResidFun,prior,priorcov,opt);
    
    % Obtain free parameter subset 
    %   - Update fixFlag to include variables with zero variance prior values
    [pinitFree,priorFree,priorcovFree] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag);

    % Test if any free parameters remain
    if(all(fixFlag==1))
        pfit = pinit;
        return;
    end

    % Define wrapper function for nLogP to use only free parameters
    nLogPFreeFun = @(pFree)(nLogPFun(getAllParams(pFree,pinit,fixFlag)));
    pfitFree = fitParams(pinitFree,nLogPFreeFun,opt);

    pfit = pinit;
    pfit(~fixFlag) = pfitFree;
end
