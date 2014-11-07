% fitErrModResid - fit error bar adjustment params given current model residuals
function [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrMod,linTransM, yresid,dydxmod,xerr,measGrpID,opt)
    % May need to add eos fitting specific options
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    % If unspecified, used identity transformation with zero offset
    if(isempty(linTransM))
        linTransM = [eye(numel(priorErrMod)) zeros(numel(priorErrMod),1)];
    end
    % If unspecified, measGrpID = '1' for all datapoints
    if(isempty(measGrpID))
        measGrpID = '1';
    end
    % If only one group specified, measGrpID is the same for all datapoints
    if(isstr(measGrpID))
        Ndat = length(yresid);
        measGrpID = cellstr(repmat(measGrpID,Ndat,1));
    end

    checkInput(pinitErrMod,priorErrMod,priorcovErrMod,linTransM,yresid,dydxmod,xerr,measGrpID);

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
        yresid,dydxmod,xerr,measGrpID)
    dimx = size(dydxmod,2);
    NpErrMod = length(pinitErrMod);
    Ndat = size(yresid,1);
    uniqID = unique(measGrpID);


    assert(length(priorErrMod)==NpErrMod,...
        'Number of error mod params must be equal in prior and init');
    assert(all(NpErrMod*[1 1]==size(priorcovErrMod)),...
        'priorcov matrix must be square with dimensions of param num');

    assert(size(yresid,2)==1,'yresid must be vertical array')

    assert(size(dydxmod,1)==Ndat,'dydxmod must have one row per datum.')
    assert(size(xerr,1)==Ndat,'xerr must have one row per datum.')

    assert(size(xerr,2)==dimx,'xerr must have one col per x dimension.')

    assert(size(linTransM,1)==dimx,...
        'linTransM must have one row per x dimension.');
    assert(size(linTransM,2)==NpErrMod+1,['linTransM must have one col per '...
        'errMod parameter plus one extra for offsets.']);

    assert(iscellstr(measGrpID),...
        'measGrpID must be a cell array of strings identifying each group');
    assert(length(measGrpID)==Ndat,'measGrpID must have one entry per datum.')

    %assert(uniqID(1)==1,'measGrpID must start with 1.');
    %assert(all(diff(uniqID)==1),'measGrpID must be an array of integers');
    %assert(dimx*length(uniqID)==NpErrMod,...
    %    'Number of error mod params must be equal to the num of measGrps*dimx');
end
