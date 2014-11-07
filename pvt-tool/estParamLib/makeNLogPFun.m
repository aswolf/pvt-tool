% makeNLogPFun - return handle for -logP fun with robust options 
%
function nLogPFun = makeNLogPFun(wtResidFun,prior,priorcov,opt)
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    assert(all(diag(priorcov)>0 | isnan(diag(priorcov))),...
        ['priorcov must have positive diagonal vals. If param is fixed, set '...
        'its diagonal to NaN.']);

    % Filter out unconstrained variables from prior weighting
    %  - unconstrained variables shown by infinite variance
    %  - NaN values mark fixed params post fitting, so ignore these too
    infVar = isinf(diag(priorcov)) | isnan(diag(priorcov));
    priorcovSub = priorcov(~infVar,~infVar);
    priorHessSub = inv(priorcovSub);
    priorSub = prior(~infVar);

    % If robust, implement as Student's T distribution
    nu = opt.robustNormParam;
    priorMismatch = @(ptest)((ptest(~infVar)-priorSub)*priorHessSub*...
        (ptest(~infVar)-priorSub)');
    if(opt.robustPrior)
        priorWt = @(ptest)(0.5*(nu+1)*sum(log(1+priorMismatch(ptest)/nu)));
    else
        priorWt = @(ptest)(0.5*priorMismatch(ptest));
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
