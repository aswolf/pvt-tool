function tests = estParamLibTest
    tests = functiontests(localfunctions);
end

function testGetAllParams_subset(testCase)
    p       = [1 2 3 4 5 6];
    fixFlag = [0 0 0 1 1 0];
    pFree = p(~fixFlag);
    
    [ptest] = getAllParams(pFree,p,fixFlag);
    verifyTrue(testCase,all(ptest==p));
end
function testGetFreeParams_subset(testCase)
    pinit   = [0 1 2 3 4];
    prior   = [5 6 7 8 9];
    priorcov = 1+[0:4;5:9;10:14;15:19;20:24];
    fixFlag = [0 1 0 0 1];
    priorcovFree = [1 3 4; 11 13 14; 16 18 19];

    [pinitFreeTest,priorFreeTest,priorcovFreeTest] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag);
    verifyTrue(testCase,all(pinitFreeTest==pinit(~fixFlag)));
    verifyTrue(testCase,all(priorFreeTest==prior(~fixFlag)));
    verifyTrue(testCase,all(reshape(priorcovFree==priorcovFreeTest,[],1)));
end
function testGetFreeParams_fixFlag(testCase)
    pinit   = [0 1 2 3 4];
    prior   = [5 6 7 8 9];
    priorcov = eye(length(prior));
    fixFlag = [1 0 0 4 0];
    
    [pinitFreeTest,priorFreeTest,priorcovFreeTest] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag);
    verifyTrue(testCase,all(pinitFreeTest==pinit(~fixFlag)));
end
%function testGetFreeParams_zeroprior(testCase)
%    keyboard;
%    pinit   = [0 1 2 3 4];
%    prior   = [5 6 7 8 9];
%    priorcov = ones(length(prior));
%    priorcov(1,1) = 0;
%    priorcov(4,4) = 0;
%    fixFlag = [0 0 0 0 0];
%    
%
%    assertError(testCase,@()(getFreeParams(pinit,prior,priorcov,fixFlag)),...
%        '''''');
%
%    [pinitFreeTest,priorFreeTest,priorcovFreeTest] = ...
%        getFreeParams(pinit,prior,priorcov,fixFlag);
%    fixFlagAll = fixFlag;
%    fixFlagAll(diag(priorcov)==0) = 1;
%    verifyTrue(testCase,all(pinitFreeTest==pinit(~fixFlagAll)));
%end
function testUpdateFreeCov_subset(testCase)
    prior = [1 2 3 4 5 6];
    fixFlagSubset = [0 0 0 1 1 1];
    pinit = prior + [.1,-.2,+.3,0,0,0];
    priorcov = diag(3*ones(size(prior)));

    [pinitFree,priorFree,priorcovFree] = ...
        getFreeParams(pinit,prior,priorcov,fixFlagSubset);

    pFree = pinitFree +10;
    pcovFree = eye(length(pFree));
    pcovFree(1,3) = .5;
    pcovFree(3,1) = .5;

    [p]    = getAllParams(pFree,prior,fixFlagSubset);
    [pcov] = updateFreeCov(pcovFree,priorcov,fixFlagSubset);

    verifyTrue(testCase,all(size(p)==size(prior)));
    verifyTrue(testCase,all(size(pcov)==size(priorcov)));
    verifyTrue(testCase,all(reshape(pcov(1:3,1:3)==pcovFree,[],1)));
    verifyTrue(testCase,all(reshape(pcov(4:6,4:6) == priorcov(4:6,4:6),[],1)));
end
function testUpdateFreeCov_nofix(testCase)
    prior = [1 2 3 4 5 6];
    fixFlag = [0 0 0 0 0 0];
    pinit = prior + [.1,-.2,+.3,-.4,-.5,-.6];
    priorcov = diag(3*ones(size(prior)));

    pcovFree = eye(size(priorcov));

    [pcov] = updateFreeCov(pcovFree,priorcov,fixFlag);

    verifyTrue(testCase,all(size(pcov)==size(priorcov)));
    verifyTrue(testCase,all(reshape(pcov==pcovFree,[],1)));
end

function testGetEstParamDefaultOpt(testCase)
    optDefault = getEstParamDefaultOpt();
    verifyTrue(testCase,isstruct(optDefault));
    verifyTrue(testCase,length(fieldnames(optDefault))>0);
end
    
function testSetDefaultOpt_empty(testCase)
    opt = [];
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);
    verifyTrue(testCase,isequal(opt,optDefault));
end
function testFitParams_uncorr(testCase)
    TOL = 1e-5;
    pinit   = [1 -2 3];
    nLogPFun = @(p)(sum(p.^2));
    opt = [];
    [pfit] = fitParams(pinit,nLogPFun,opt);
    verifyTrue(testCase,all(abs(pfit-[0 0 0]) < TOL));
end
function testFitParams_corr(testCase)
    TOL = 1e-5;
    pinit   = [1 2];
    pcov = [1 .95; .95 1];
    phess = inv(pcov);
    nLogPFun = @(p)(p*phess*p');
    opt = [];
    [pfit] = fitParams(pinit,nLogPFun,opt);
    verifyTrue(testCase,all(abs(pfit-[0 0]) < TOL));
end
function testFitParams_corrNd_NOT(testCase)
    verifyTrue(testCase,false);
    % creating a random valid covariance matrix (positive semi definite)
    % is actually non trivial!

    TOL = 1e-5;
    Nd = 5;
    rndSeed = 42+1;
    rng(rndSeed);
    [Q, R] = qr(100*randn(Nd));
    Q = Q*diag(sign(diag(R)));
    if(abs(1+det(Q))<TOL)
        col = Q(:,1);
        Q(:,1) = -1*col;
    end
    rotv = Q*eye(Nd);
    rotv*eye(Nd)*inv(rotv);

    %priorSig = ones(1,Nd);
    %pcov = .9*(2*rand(Nd)-1);
    phess = randn(Nd);
    for(i=1:Nd)
        phess(i,i) = abs(phess(i,i));
    end
    pinit = randn(1,Nd);

    %NOTE: prior is necessary to prevent numerical overflow!
    phess = inv(pcov);
    nLogPFun = @(p)(p*phess*p' + sum((p./priorSig).^2));
    nLogPFun(pinit);
    opt = [];
    [pfit] = fitParams(pinit,nLogPFun,opt);
    %verifyTrue(testCase,all(abs(pfit-[0 0]) < TOL));
end

function testEstParamCov_uncorrQuad(testCase)
    TOL = 1e-5;
    nLogPFun = @(p)(0.5*sum(p.^2));
    pfit = [0 0 0 0];
    opt = [];
    fixFlag = [];
    [pfitcov] = estParamCov(nLogPFun,pfit,fixFlag,opt);

    dpfitcov = pfitcov - eye(length(pfit));
    verifyTrue(testCase,all(abs(dpfitcov(:)) < TOL));
end
function testEstParamCov_uncorrSclQuad(testCase)
    TOL = 1e-5;
    pSig = [1 3 10 30];
    nLogPFun = @(p)(0.5*sum((p./pSig).^2));
    pfit = [0 0 0 0];
    opt = [];
    fixFlag = [];
    [pfitcov] = estParamCov(nLogPFun,pfit,fixFlag,opt);

    dpfitcov = pfitcov - diag(pSig.^2);
    verifyTrue(testCase,all(abs(dpfitcov(:)) < TOL));
end

function testEstParamCov_corr2d(testCase)
    TOL = 1e-5;
    pcorr = [1 .9; .9 1];
    pcov = pcorr;
    phess = inv(pcov);
    nLogPFun = @(p)(0.5*p*phess*p');
    pfit = [0 0];
    opt = [];
    fixFlag = [];
    [pfitcov] = estParamCov(nLogPFun,pfit,fixFlag,opt)

    dpfitcov = pfitcov - pcov;
    verifyTrue(testCase,all(abs(dpfitcov(:)) < TOL));
end
function testEstParamCov_corrNd_NOT(testCase)
    verifyTrue(testCase,false);
end

function testFitFreeParamsWithPrior_NOT(testCase)
    verifyTrue(testCase,false);
end
function testMkNLogPFun_optDefaults(testCase)
    optIn.NfitIter = 5;
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(optIn,optDefault);
    optNameList = fieldnames(optDefault);
    optInList = fieldnames(opt);
    otherFields = setdiff(optNameList,optInList);
    verifyTrue(testCase,isempty(otherFields));
end
function testMkNLogPFun_dist_NOT(testCase)
    % Test that distribution of costfunction values matches a chi sqr distbn
    verifyTrue(testCase,false);
end
function testMkNLogPFun_distRobust_NOT(testCase)
    % Test that distribution of costfunction values matches a multivariate
    % students t cdf
    verifyTrue(testCase,false);
end

function testCovToCorr_simple(testCase)
    TOL = 1e-6;
    pcorr0 = [1 .5; .5 1];
    perr0  = [3 10];
    perrM0 = perr0(:)*perr0(:)';
    pcov = perrM0.*pcorr0;
    
    [perr,pcorr] = covToCorr(pcov);

    verifyTrue(testCase, all(perr(:)-perr0(:)<TOL));
    verifyTrue(testCase, all(pcorr(:)-pcorr0(:)<TOL));
end
function testCovToCorr_InfElem(testCase)
    TOL = 1e-6;
    pcorr0 = [1 .5 0 ; .5 1 0; 0 0 1];
    perr0  = [3 10 Inf];
    perrM0 = perr0(:)*perr0(:)';
    pcov = perrM0.*pcorr0;
    pcov(isnan(pcov)) = 0;
    
    [perr,pcorr] = covToCorr(pcov);

    verifyTrue(testCase, all(perr(:)==perr0(:)));
    verifyTrue(testCase, all(pcorr(:)-pcorr0(:)<TOL));
end
function testCorrToCov_simple(testCase)
    TOL = 1e-6;
    pcorr0 = [1 .5; .5 1];
    perr0  = [3 10];
    perrM0 = perr0(:)*perr0(:)';
    pcov0 = perrM0.*pcorr0;
    
    pcov = corrToCov(perr0,pcorr0);

    verifyTrue(testCase, all(pcov(:)-pcov0(:)<TOL));
end
function testCorrToCov_InfElem(testCase)
    pcorr0 = [1 .5 0 ; .5 1 0; 0 0 1];
    perr0  = [3 10 Inf];
    perrM0 = perr0(:)*perr0(:)';
    pcov0 = perrM0.*pcorr0;
    pcov0(isnan(pcov0)) = 0;
    
    pcov = corrToCov(perr0,pcorr0);

    verifyTrue(testCase, all(pcov(:)==pcov0(:)));
end
function testCorrToCov_InfElemWithCorr(testCase)
    pcorr0 = [1 .5 .1 ; .5 1 0; .1 0 1];
    perr0  = [3 10 Inf];
    perrM0 = perr0(:)*perr0(:)';
    pcov0 = perrM0.*pcorr0;
    pcov0(isnan(pcov0)) = 0;
    pcov0(isinf(pcov0)) = 0;
    infInd = find(isinf(perr0));
    pcov0(sub2ind(size(pcov0),infInd,infInd)) = Inf;

    
    pcov = corrToCov(perr0,pcorr0);

    verifyTrue(testCase, all(pcov(:)==pcov0(:)));
end
function testFitErrModResid_1dPolyDist(testCase)
    seed=42;
    rng(seed)

    TOLWID = 0.3;
    Ndraw = 100;

    %NOTE: errorbar MUST be small enough that linear approximation is
    %reasonable
    xerrMag = .05;
    ptrueErrMod = [-0.4];
    priorcovErrMod = [1^2];

    relDev = zeros(Ndraw,1);
    for(i=1:Ndraw)
        % Test errMod for polynomial
        pmod = [0.2 1.5 0];
        dpmod = polyder(pmod);

        xmod = linspace(-3,3,30)';
        ymod = polyval(pmod,xmod);

        xerr = xerrMag*ones(size(xmod));
        xdev = exp(ptrueErrMod)*xerr.*randn(size(xerr));
        xobs = xmod + xdev;

        %plot(xmod,ymod,'ro',xobs,ymod,'kx')
        yresid = ymod - polyval(pmod,xobs);
        dydxobs = polyval(dpmod,xobs);

        measGrpID = cellstr(num2str(ones(size(yresid))));

        pinitErrMod = [0];
        priorErrMod = [0];
        opt = [];
        linTransM = [];

        [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
            priorErrMod,priorcovErrMod,linTransM,yresid,dydxobs,xerr,measGrpID,opt);

        [pfitcovErrMod] = estParamCov(nLogPFun,pfitErrMod,[],opt);

        irelDev = (ptrueErrMod-pfitErrMod)/sqrt(pfitcovErrMod);
        relDev(i) = irelDev;
    end

    %hist(relDev,10)
    %mean(relDev)
    %std(relDev)

    verifyTrue(testCase,abs(mean(relDev))/std(relDev) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    verifyTrue(testCase,abs(log(std(relDev))) < TOLWID ,...
        ['Distribution of normalized errMod params must have width close to '...
        'truth (ie. within TOLWID of 1)']);

end
function testFitErrModResid_1dPolyMultiMeasGrp_NOT(testCase)
    verifyTrue(testCase,false ,'Multiple measGrpIDs not yet implemented');
end

function testFitErrModResid_2dPolyFunScl(testCase)
    seed=2*42;
    rng(seed);

    TOLWID = 0.1;
    Ndraw = 100;

    x1errMag = 0.03;
    x2errMag = 0.01;
    ptrueErrModScl = -1;
    ptrueErrMod = ptrueErrModScl*[1 1];
    pinitErrMod = [0 0];
    priorErrMod = [0 0];
    priorcovErrMod = diag([3^2 3^2]);

    opt = [];
    linTransM = [];

    relDev = zeros(Ndraw,1);

    % Test errMod for polynomial
    %             1      2     3   4   5   6
    %           x1^2    x2^2 x1*x2 x1  x2  1
    pmod     = [0.2     -.2    0  -1.5 1.5 0];
    %dydx1 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + 2*pmod(1)*x1 + x2*pmod(3) + pmod(4)
    %dydx2 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + x1*pmod(3)   + 2*pmod(2)*x2 +  pmod(5)
    dpmoddx1 = [0 0 0 2*pmod(1)   pmod(3) pmod(4)];
    dpmoddx2 = [0 0 0   pmod(3) 2*pmod(2) pmod(5)];
    
    x1vec = linspace(-3  ,3  ,10);
    x2vec = linspace(-2.5,3.5,10);
    [x1M,x2M] = ndgrid(x1vec,x2vec);
    x12 = [x1M(:) x2M(:)];
    Ndat = size(x12,1);

    for(i=1:Ndraw)
        xerr = ones(size(x12,1),1)*[x1errMag x2errMag];
        xdev = repmat(exp(ptrueErrMod),Ndat,1).*xerr.*randn(size(xerr));
        xobs = x12 + xdev;

        designM = [x12(:,1).^2 x12(:,2).^2 x12(:,1).*x12(:,2) ...
            x12(:,1) x12(:,2) ones(size(x12,1),1)];

        designObsM = [xobs(:,1).^2 xobs(:,2).^2 xobs(:,1).*xobs(:,2) ...
            xobs(:,1) xobs(:,2) ones(size(xobs,1),1)];

        ymod    = designM*pmod';

        yresid  = ymod-designObsM*pmod';
        dydx1obs= designObsM*dpmoddx1';
        dydx2obs= designObsM*dpmoddx2';
        dydxobs = [dydx1obs(:) dydx2obs(:)];

        yM = reshape(ymod,size(x1M));

        measGrpID = cellstr(num2str(ones(size(yresid))));


        [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
            priorErrMod,priorcovErrMod,linTransM,yresid,dydxobs,xerr,measGrpID,opt);

        nLogPFunScl=@(scl)(nLogPFun(scl*[1 1]));
        pfitErrModScl = fminunc(nLogPFunScl,0,optimset('Display','off','LargeScale','off'));
        [pfitcovErrMod] = estParamCov(nLogPFunScl,pfitErrModScl,[],opt);

        irelDev = (ptrueErrModScl-pfitErrModScl)/sqrt(pfitcovErrMod);
        relDev(i) = irelDev;

        %irelDev = (ptrueErrMod-pfitErrMod)./sqrt(diag(pfitcovErrMod)');
        %relDev(i,:) = irelDev;
    end
    %keyboard;
    %hist(relDev,10)
    %mean(relDev)
    %std(relDev)

    assert(abs(mean(relDev))/std(relDev) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    assert(abs(log(std(relDev))) < TOLWID ,...
        ['Distribution of normalized errMod params must have width close to '...
        'truth (ie. within TOLWID of 1)']);

end
function testFitErrModResid_2dPolyLinTransEquiv(testCase)
    TOL = 1e-4;
    seed=2*42;
    rng(seed);

    x1errMag = 0.03;
    x2errMag = 0.01;
    ptrueErrModScl = -1;
    ptrueErrMod = ptrueErrModScl*[1 1];
    pinitErrMod = [0 0];
    priorErrMod = [0 0];
    priorcovErrMod = diag([1e3^2 1e3^2]);

    opt = [];

    % Test errMod for polynomial
    %             1      2     3   4   5   6
    %           x1^2    x2^2 x1*x2 x1  x2  1
    pmod     = [0.2     -.2    0  -1.5 1.5 0];
    %dydx1 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + 2*pmod(1)*x1 + x2*pmod(3) + pmod(4)
    %dydx2 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + x1*pmod(3)   + 2*pmod(2)*x2 +  pmod(5)
    dpmoddx1 = [0 0 0 2*pmod(1)   pmod(3) pmod(4)];
    dpmoddx2 = [0 0 0   pmod(3) 2*pmod(2) pmod(5)];
    
    x1vec = linspace(-3  ,3  ,10);
    x2vec = linspace(-2.5,3.5,10);
    [x1M,x2M] = ndgrid(x1vec,x2vec);
    x12 = [x1M(:) x2M(:)];
    Ndat = size(x12,1);

    xerr = ones(size(x12,1),1)*[x1errMag x2errMag];
    xdev = repmat(exp(ptrueErrMod),Ndat,1).*xerr.*randn(size(xerr));
    xobs = x12 + xdev;

    designM = [x12(:,1).^2 x12(:,2).^2 x12(:,1).*x12(:,2) ...
        x12(:,1) x12(:,2) ones(size(x12,1),1)];

    designObsM = [xobs(:,1).^2 xobs(:,2).^2 xobs(:,1).*xobs(:,2) ...
        xobs(:,1) xobs(:,2) ones(size(xobs,1),1)];

    ymod    = designM*pmod';

    yresid  = ymod-designObsM*pmod';
    dydx1obs= designObsM*dpmoddx1';
    dydx2obs= designObsM*dpmoddx2';
    dydxobs = [dydx1obs(:) dydx2obs(:)];

    yM = reshape(ymod,size(x1M));

    measGrpID = cellstr(num2str(ones(size(yresid))));


    [tr nLogPFun] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrMod,[],yresid,dydxobs,xerr,measGrpID,opt);

    nLogPFunScl=@(scl)(nLogPFun(scl*[1 1]));

    %compare with linTransM
    linTransM = [1 0.5 0; 1 -0.5 0];
    priorcovErrModT = diag([1e3^2 1e-3^2]);
    [tr nLogPFunT] = fitErrModResid(pinitErrMod,...
        priorErrMod,priorcovErrModT,linTransM,yresid,dydxobs,xerr,measGrpID,opt);

    sclMag = linspace(-2,0,30);
    nLogPScl  = zeros(size(sclMag));
    nLogPTrans= zeros(size(sclMag));
    for(i=1:length(sclMag))
        nLogPScl(i) = nLogPFunScl(sclMag(i));
        nLogPTrans(i) = nLogPFunT([sclMag(i)  0]);
    end
    %plot(sclMag,nLogPScl,'k-',sclMag,nLogPTrans,'ro')

    verifyTrue(testCase,all(abs(nLogPScl-nLogPTrans)<TOL),...
        'Linear transformation must produce correct result when x1=x2');

end
function testFitErrModResid_2dPolyLinTransFitScl(testCase)
    seed=2*42;
    rng(seed);


    TOLWID = 0.3;
    Ndraw = 30;

    x1errMag = 0.003;
    x2errMag = 0.001;

    linTransM = [1 0.5 0; 1 -0.5 0];
    ptrueErrMod = [+.2 +.04];
    ptrueErrMod = [+.2 +0];

    pinitErrMod = [0 0];
    priorErrMod = [0 0];
    priorcovErrMod = diag([3^2 .001^2]);
    %priorcovErrMod = diag([3^2 3^2]);


    relDev = zeros(Ndraw,2);
    absDev = zeros(Ndraw,2);
    fitList = zeros(Ndraw,2);

    % Test errMod for polynomial
    %             1      2     3   4   5   6
    %           x1^2    x2^2 x1*x2 x1  x2  1
    pmod     = [0.2     -.3    0  -1.5 -3 0];
    %dydx1 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + 2*pmod(1)*x1 + x2*pmod(3) + pmod(4)
    %dydx2 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + x1*pmod(3)   + 2*pmod(2)*x2 +  pmod(5)
    dpmoddx1 = [0 0 0 2*pmod(1)   pmod(3) pmod(4)];
    dpmoddx2 = [0 0 0   pmod(3) 2*pmod(2) pmod(5)];
    
    x1vec = linspace(-3  ,3  ,10);
    x2vec = linspace(-2.5,3.5,10);
    [x1M,x2M] = ndgrid(x1vec,x2vec);
    x12 = [x1M(:) x2M(:)];
    Ndat = size(x12,1);
    xerrFac = exp(linTransM*[ptrueErrMod 1]')';

    for(i=1:Ndraw)
        xerr = ones(size(x12,1),1)*[x1errMag x2errMag];
        xdev = repmat(xerrFac,Ndat,1).*xerr.*randn(size(xerr));
        xobs = x12 + xdev;

        designM = [x12(:,1).^2 x12(:,2).^2 x12(:,1).*x12(:,2) ...
            x12(:,1) x12(:,2) ones(size(x12,1),1)];

        designObsM = [xobs(:,1).^2 xobs(:,2).^2 xobs(:,1).*xobs(:,2) ...
            xobs(:,1) xobs(:,2) ones(size(xobs,1),1)];

        ymod    = designM*pmod';

        yresid  = ymod-designObsM*pmod';
        dydx1obs= designObsM*dpmoddx1';
        dydx2obs= designObsM*dpmoddx2';
        dydxobs = [dydx1obs(:) dydx2obs(:)];

        yM = reshape(ymod,size(x1M));

        measGrpID = cellstr(num2str(ones(size(yresid))));

        opt = [];

        [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
            priorErrMod,priorcovErrMod,linTransM,yresid,dydxobs,xerr,...
            measGrpID,opt);

        %sclMag = linspace(-.6,.6,30);
        %nLogPScl  = zeros(size(sclMag));
        %nLogPTrans= zeros(size(sclMag));
        %for(i=1:length(sclMag))
        %    nLogPTrans(i) = nLogPFun([sclMag(i)  0]);
        %end
        %plot(sclMag,nLogPTrans,'k-')

        [pfitcovErrMod] = estParamCov(nLogPFun,pfitErrMod,[],opt);
        iabsDev = (ptrueErrMod-pfitErrMod);
        irelDev = (ptrueErrMod-pfitErrMod)./sqrt(diag(pfitcovErrMod))';

        absDev(i,:) = iabsDev;
        relDev(i,:) = irelDev;
        fitList(i,:) = pfitErrMod;
    end

    verifyTrue(testCase,abs(mean(relDev(:,1)))/std(relDev(:,1)) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    verifyTrue(testCase,abs(log(std(relDev(:,1)))) < TOLWID ,...
        ['Distribution of normalized errMod param for avg magnitude must '...
        'have width close to truth (ie. within TOLWID of 1)']);
end
function testFitErrModResid_2dPolyLinTransFit_NOT(testCase)
    seed=2*42;
    rng(seed);

    verifyTrue(testCase,false,['2dPolyLinTransFit not yet implemented.' ...
        'Cant get diff parameter yet']);

    TOLWID = 0.3;
    Ndraw = 300;

    x1errMag = 0.0003;
    x2errMag = 0.0001;

    linTransM = [1 0.5 0; 1 -0.5 0];
    ptrueErrMod = [+.2 -.04];

    pinitErrMod = [0 0];
    pinitErrMod = ptrueErrMod;
    priorErrMod = [0 0];
    priorcovErrMod = diag([.3^2 .1^2]);
    %priorcovErrMod = diag([3^2 3^2]);


    relDev = zeros(Ndraw,2);
    absDev = zeros(Ndraw,2);
    fitList = zeros(Ndraw,2);

    % Test errMod for polynomial
    %             1      2     3   4   5   6
    %           x1^2    x2^2 x1*x2 x1  x2  1
    pmod     = [0.2     -.3    .1  -1.5 -3 0];
    %dydx1 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + 2*pmod(1)*x1 + x2*pmod(3) + pmod(4)
    %dydx2 = 0*x1^2 + 0*x2^2 + 0*x1*x2 + x1*pmod(3)   + 2*pmod(2)*x2 +  pmod(5)
    dpmoddx1 = [  0      0     0 2*pmod(1)   pmod(3) pmod(4)];
    dpmoddx2 = [  0      0     0   pmod(3) 2*pmod(2) pmod(5)];
    
    x1vec = linspace(-3  ,3  ,10);
    x2vec = linspace(-2.5,3.5,10);
    [x1M,x2M] = ndgrid(x1vec,x2vec);
    x12 = [x1M(:) x2M(:)];
    Ndat = size(x12,1);
    xerrFac = exp(linTransM*ptrueErrMod')';

    for(i=1:Ndraw)
        xerr = ones(size(x12,1),1)*[x1errMag x2errMag];
        xdev = repmat(xerrFac,Ndat,1).*xerr.*randn(size(xerr));
        xobs = x12 + xdev;

        designM = [x12(:,1).^2 x12(:,2).^2 x12(:,1).*x12(:,2) ...
            x12(:,1) x12(:,2) ones(size(x12,1),1)];

        designObsM = [xobs(:,1).^2 xobs(:,2).^2 xobs(:,1).*xobs(:,2) ...
            xobs(:,1) xobs(:,2) ones(size(xobs,1),1)];

        ymod    = designM*pmod';

        yresid  = ymod-designObsM*pmod';
        dydx1obs= designObsM*dpmoddx1';
        dydx2obs= designObsM*dpmoddx2';
        dydxobs = [dydx1obs(:) dydx2obs(:)];

        yM = reshape(ymod,size(x1M));

        measGrpID = cellstr(num2str(ones(size(yresid))));

        opt = [];

        [pfitErrMod nLogPFun] = fitErrModResid(pinitErrMod,...
            priorErrMod,priorcovErrMod,linTransM,yresid,dydxobs,xerr,...
            measGrpID,opt);

        %x1Sv = linspace(-.5,.5,21);
        %x2Sv = linspace(-.1,.1,23);
        %[x1SM,x2SM] = meshgrid(x1v,x2v);
        %nLogPSM = zeros(length(x1Sv),length(x2Sv));
        %for(i=1:length(x1Sv))
        %    for(j=1:length(x2Sv))
        %        nLogPSM(i,j) = nLogPFun([x1Sv(i),x2Sv(j)]);
        %    end
        %end
        %contour(x2Sv,x1Sv,nLogPSM-nLogPFun(pfitErrMod),[ 0 1 2])
        %hold on;
        %plot(pfitErrMod(2),pfitErrMod(1),'kx',ptrueErrMod(2),ptrueErrMod(1),'ro')
        %hold off;

        %sclMag = linspace(-.6,.6,30);
        %nLogPScl  = zeros(size(sclMag));
        %nLogPTrans= zeros(size(sclMag));
        %for(i=1:length(sclMag))
        %    nLogPTrans(i) = nLogPFun([sclMag(i)  0]);
        %end
        %plot(sclMag,nLogPTrans,'k-')

        [pfitcovErrMod] = estParamCov(nLogPFun,pfitErrMod,[],opt);
        iabsDev = (ptrueErrMod-pfitErrMod);
        irelDev = (ptrueErrMod-pfitErrMod)./sqrt(diag(pfitcovErrMod))';

        absDev(i,:) = iabsDev;
        relDev(i,:) = irelDev;
        fitList(i,:) = pfitErrMod;
    end

    %ptrueErrMod
    %hist(fitList(:,1),10)
    %hist(fitList(:,2),10)

    %hist(relDev(:,1),10)
    %hist(relDev(:,2),10)

    verifyTrue(testCase,abs(mean(relDev(:,1)))/std(relDev(:,1)) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    verifyTrue(testCase,abs(log(std(relDev(:,1)))) < TOLWID ,...
        ['Distribution of normalized errMod param for avg magnitude must '...
        'have width close to truth (ie. within TOLWID of 1)']);
end
