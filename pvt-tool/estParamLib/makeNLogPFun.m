% makeNLogPFun - return handle for -logP fun with robust options 
%
function nLogPFun = makeNLogPFun(wtResidFun,prior,priorcov,opt)
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    % Filter out unconstrained variables from prior weighting
    %  - unconstrained variables shown by infinite variance
    infVar = isinf(diag(priorcov));
    priorcovSub = priorcov(~infVar,~infVar);
    priorHessSub = inv(priorcovSub);
    priorSub = prior(~infVar);
    ptestSub = ptest(~infVar);

    % If robust, implement as Student's T distribution
    nu = opt.robustNormParam;
    priorMismatch = (ptestSub-priorSub)*priorHessSub*(ptestSub-priorSub)';
    if(opt.robustPrior)
        priorWt = @(ptest)(0.5*(nu+1)*sum(log(1+priorMismatch/nu)));
    else
        priorWt = @(ptest)(0.5*priorMismatch);
    end

    % Construct total posterior Function nLogPFun
    %   -depends on whether it is robust or not
    if(opt.robustFit)
        nLogPFun = @(ptest)(0.5*(nu+1)*sum(log(1+wtResidFun(ptest).^2/nu)) + ...
            priorWt(ptest));
    else
        nLogPFun = @(ptest)(0.5*sum(wtResidFun(ptest).^2) + ...
            priorWt(ptest));
    end
end
