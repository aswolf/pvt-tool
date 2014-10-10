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
function testFitParams_corrNd(testCase)
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
function testEstParamCov_corrNd(testCase)
    verifyTrue(testCase,false);
end

function testFitFreeParamsWithPrior(testCase)
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
function testMkNLogPFun_dist(testCase)
    % Test that distribution of costfunction values matches a chi sqr distbn
    verifyTrue(testCase,false);
end
function testMkNLogPFun_distRobust(testCase)
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
