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
function testGetFreeParams_zeroprior(testCase)
    pinit   = [0 1 2 3 4];
    prior   = [5 6 7 8 9];
    priorcov = ones(length(prior));
    priorcov(1,1) = 0;
    priorcov(4,4) = 0;
    fixFlag = [0 0 0 0 0];
    
    [pinitFreeTest,priorFreeTest,priorcovFreeTest] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag);
    fixFlagAll = fixFlag;
    fixFlagAll(diag(priorcov)==0) = 1;
    verifyTrue(testCase,all(pinitFreeTest==pinit(~fixFlagAll)));
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

function testEstParamCov_uncorr(testCase)
    TOL = 1e-5;
    nLogPFun = @(p)(0.5*sum(p.^2));
    pfit = [0 0 0 0];
    opt = [];
    [pfitcov] = estParamCov(nLogPFun,pfit,opt);

    dpfitcov = pfitcov - eye(length(pfit));
    verifyTrue(testCase,all(abs(dpfitcov(:)) < TOL));
end
function testEstParamCov_uncorrScl(testCase)
    TOL = 1e-5;
    pSig = [1 3 10 30];
    nLogPFun = @(p)(0.5*sum((p./pSig).^2));
    pfit = [0 0 0 0];
    opt = [];
    [pfitcov] = estParamCov(nLogPFun,pfit,opt);

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
    [pfitcov] = estParamCov(nLogPFun,pfit,opt)

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
