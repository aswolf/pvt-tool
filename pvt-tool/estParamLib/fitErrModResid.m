% fitErrModResid - fit error bar adjustment params given current model residuals
function [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrMod,linTransM, yresid,dydxmod,xerr,opt)
    priorcovErrMod = squeeze(priorcovErrMod);
    % May need to add eos fitting specific options
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    % If unspecified, used identity transformation with zero offset
    if(isempty(linTransM))
        linTransM = [eye(numel(priorErrMod)) zeros(numel(priorErrMod),1)];
    end

    checkInput(pinitErrMod,priorErrMod,priorcovErrMod,linTransM,yresid,dydxmod,xerr);

    pinitErrMod = pinitErrMod(:)';
    Ndat = length(yresid);
    
    sqrErrMod = @(perrmod)(sum((dydxmod.*xerr.*(ones(Ndat,1)*exp((linTransM*[perrmod 1]')'))).^2,2));

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

function checkInput(pinitErrMod,priorErrMod,priorcovErrMod,linTransM,...
        yresid,dydxmod,xerr)
    dimx = size(dydxmod,2);
    Ndat = size(yresid,1);
    NpErrMod = length(pinitErrMod);

    assert(all(length(priorErrMod)==length(pinitErrMod)),...
        'Number of error mod params must be equal in prior and init');
    assert(all(size(priorcovErrMod)==NpErrMod*[1 1]),...
        'priorcov matrix must provide square with dimensions of param num');

    assert(size(yresid,2)==1,'yresid must be vertical array')

    assert(size(dydxmod,1)==Ndat,'dydxmod must have one row per datum.')
    assert(size(xerr,1)==Ndat,'xerr must have one row per datum.')

    assert(size(xerr,2)==dimx,'xerr must have one col per x dimension.')

    assert(size(linTransM,1)==dimx,...
        'linTransM must have one row per x dimension.');
    assert(size(linTransM,2)==NpErrMod+1,['linTransM must have one col per '...
        'errMod parameter plus one extra for offsets.']);
end
