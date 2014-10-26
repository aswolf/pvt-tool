function [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrMod,yresid,dydxmod,xerr,measGrpID,opt)
    % May need to add eos fitting specific options
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    checkInput(pinitErrMod,priorErrMod,priorcovErrMod,yresid,dydxmod,xerr,measGrpID);

    pinitErrMod = pinitErrMod(:)';
    Ndat = length(yresid);

    
    sqrErrMod = @(perrmod)(sum((dydxmod.*xerr.*(ones(Ndat,1)*exp(perrmod))).^2,2));

    priorHess = inv(priorcovErrMod);
    % If robust, implement as Student's T distribution
    nu = opt.robustNormParam;
    priorMismatch = @(ptest)((ptest-priorErrMod)*priorHess*...
        (ptest-priorErrMod)');
    if(opt.robustPrior)
        priorWt = @(ptest)(0.5*(nu+1)*sum(log(1+priorMismatch(ptest)/nu)));
    else
        priorWt = @(ptest)(0.5*priorMismatch(ptest));
    end
    
    nLogPFun = @(perrmod)(0.5*sum(yresid.^2./sqrErrMod(perrmod) + ...
        log(sqrErrMod(perrmod)))+priorWt(perrmod));

    [pfitErrMod] = fitParams(pinitErrMod,nLogPFun,opt);
end

function checkInput(pinitErrMod,priorErrMod,priorcovErrMod,...
        yresid,dydxmod,xerr,measGrpID)
    uniqID = unique(measGrpID);
    assert(size(yresid,2)==1,'yresid must be vertical array')
    assert(size(dydxmod,1)==size(yresid,1),...
        'dydxmod must have one row per datum.')
    assert(size(xerr,1)==size(yresid,1),'xerr must have one row per datum.')
    assert(size(measGrpID,1)==size(yresid,1),'measGrpID must have same len as yresid.')
    assert(size(xerr,2)==size(dydxmod,2),...
        'xerr and dydxmod must have one col per independent var.')
    assert(uniqID(1)==1,'measGrpID must start with 1.');
    assert(all(diff(uniqID)==1),'measGrpID must be an array of integers');
    assert(size(dydxmod,2)*length(uniqID)==length(pinitErrMod),...
        'Number of error mod params must be equal to the number of runs*3');
    assert(length(priorErrMod)==length(pinitErrMod),...
        'Number of error mod params must be equal in prior and init');
    assert(all(length(priorErrMod)*[1 1]==size(priorcovErrMod)),...
        'priorcov matrix must be square with dimensions of param num');
end
