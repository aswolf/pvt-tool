function tests = fitLibTest
    tests = functiontests(localfunctions);
end

function testFitColdCompressData_zeroMeasErr(testCase)
    TOL = 1e-4;

    pEos = [150 250 4.5];
    coldEosFun = @VinetEos;
    priorEos = pEos;
    priorcovEos = diag(Inf*[1 1 1]);

    V = pEos(1)*linspace(.7,1,10)';
    Pmod = coldEosFun(V,pEos);
    PerrTot = 1e-6*ones(size(V));
    P = Pmod;
    pInitEos = pEos + [10 -15 -.5];
    %PInit = coldEosFun(V,pInitEos);
    fixFlag = [];
    opt = [];
    [pfit pfitcov] = fitColdCompressData(pInitEos,fixFlag,...
        priorEos,priorcovEos,coldEosFun,P,V,PerrTot,opt);
    relErr = (pfit-pEos)./pEos;

    verifyTrue(testCase,all(abs(relErr)<TOL));
end
function testFitColdCompressData_fixSubset(testCase)
    TOL = 1e-4;

    pEos = [150 250 4.5];
    coldEosFun = @VinetEos;
    priorEos = pEos;
    priorcovEos = diag(Inf*[1 1 1]);

    V = pEos(1)*linspace(.7,1,10)';
    Pmod = coldEosFun(V,pEos);
    PerrTot = 1e-6*ones(size(V));
    P = Pmod;
    pInitEos = pEos + [0 -15 -.5];
    %PInit = coldEosFun(V,pInitEos);
    fixFlag = [0 0 0];
    fixFlag(1) = 1;
    opt = [];
    [pfit pfitcov] = fitColdCompressData(pInitEos,fixFlag,...
        priorEos,priorcovEos,coldEosFun,P,V,PerrTot,opt);
    relErr = (pfit-pEos)./pEos;

    verifyTrue(testCase,all(abs(relErr)<TOL));
end
function testFitColdCompressData_credRegion(testCase)

    rndSeed = 42;
    rng(rndSeed)
    Perr = 0.5;

    pEos = [150 250 4.5];
    coldEosFun = @VinetEos;
    priorEos = pEos;
    priorcovEos = diag(Inf*[1 1 1]);

    V = pEos(1)*linspace(.7,1,10)';
    Pmod = coldEosFun(V,pEos);
    PerrTot = Perr*ones(size(V));
    dof = length(pEos);
    %probLvl = diff(normcdf([-1 1]));
    probLvl = 0.5;
    chi2Thresh = invchi2(probLvl,dof);
    pInitEos = pEos;
    fixFlag = [];
    opt.NfitIter  = 1;


    Niter = 30;
    costFunHist = zeros(Niter,1);
    chi2MinHist = zeros(Niter,1);
    chi2TrueHist = zeros(Niter,1);
    pEosHist = zeros(Niter,length(pEos));
    pcovHist = zeros(Niter,length(pEos),length(pEos));
    for(i=1:Niter)
        iP = Pmod+PerrTot.*randn(size(V));
        %PInit = coldEosFun(V,pInitEos);
        [ipfit ipfitcov inLogPFun] = fitColdCompressData(pInitEos,fixFlag,...
            priorEos,priorcovEos,coldEosFun,iP,V,PerrTot,opt);
        costFunHist(i) = (ipfit-pEos)*inv(ipfitcov)*(ipfit-pEos)';
        chi2MinHist(i) = 2*inLogPFun(ipfit);
        chi2TrueHist(i) = 2*inLogPFun(pEos);
        pEosHist(i,:) = ipfit;
        pcovHist(i,:,:) = ipfitcov;
    end

    costFunHists = sort(costFunHist);
    chi2MinHists = sort(chi2MinHist);
    chi2TrueHists = sort(chi2TrueHist);

    indMed = round(.5*Niter);

    Ntry = 3000;
    chi2draw = chi2rnd(dof,Niter,Ntry);
    chi2draws = sort(chi2draw);
    chi2Thresh = quantile(chi2draws(indMed,:),0.95);
    popFrac = mean(costFunHist<chi2Thresh);

    verifyTrue(testCase,popFrac > 0.5);
end

function testFitHotCompressData_zeroMeasErr(testCase)
    TOL = 1e-3;

    T0 = 300;
    V0 = 162.373;
    K0 = 258.4;
    KP0= 4.10;

    Natom = 4*5;
    Tdeb0 = 940;
    gam0  = 1.55;
    q = 1.1;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 q 1.0];

    VCold = V0*linspace(.7,1,9)';
    TCold = 300*ones(size(VCold));
    THot = [...
        [1000:1500:5500]'; ...
        [1000:1500:5500]'; ...
        ];
    VHot  = [...
        VCold(1)*ones(4,1); ...
        VCold(9)*ones(4,1); ...
        ];
    V = [VCold;VHot];
    T = [TCold;THot];


    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyePowerLaw;
    hotExtraInputs = {Natom, debyeDerivsFun};
    hotEosFun  = @(V,T,V0,T0,pHotEos,hotExtraInputs)...
        (MieGrunDebyeHotEos(V,T,V0,T0,pHotEos,hotExtraInputs{:}));

    addedThermPressFun = [];
    [Pmod,KTmod,Cvmod,gammod] = calcPressThermAddEos(V(:),T(:),T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);

    %scatter(Pmod,V,50,T,'o')

    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    priorEos = pEos;
    priorcovEos = diag(Inf*ones(size(pEos)));

    pInitEos = pEos + [10 -15 -.5 0 -.2 .3 0];
    fixFlag = zeros(size(priorEos));
    fixFlag(4)   = 1;
    fixFlag(end) = 1;
    %priorcovEos(end,end) = 0;

    PerrTot = 1e-6*ones(size(V));
    P = Pmod;

    %PInit = coldEosFun(V,pInitEos);
    opt = [];
    opt.NfitIter = 1;
    [pfit pfitcov] = fitHotCompressData(pInitEos,fixFlag,T0,...
        NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun,P,V,T,PerrTot,opt);
    relErr = (pfit-pEos)./pEos;

    verifyTrue(testCase,all(abs(relErr)<TOL));
end
function testFitHotCompressData_fixSubset(testCase)
    TOL = 1e-3;

    T0 = 300;
    V0 = 162.373;
    K0 = 258.4;
    KP0= 4.10;

    Natom = 4*5;
    Tdeb0 = 940;
    gam0  = 1.55;
    q = 1.1;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 q 1.0];

    VCold = V0*linspace(.7,1,9)';
    TCold = 300*ones(size(VCold));
    THot = [...
        [1000:1500:5500]'; ...
        [1000:1500:5500]'; ...
        ];
    VHot  = [...
        VCold(1)*ones(4,1); ...
        VCold(9)*ones(4,1); ...
        ];
    V = [VCold;VHot];
    T = [TCold;THot];


    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyePowerLaw;
    hotExtraInputs = {Natom, debyeDerivsFun};
    hotEosFun  = @(V,T,V0,T0,pHotEos,hotExtraInputs)...
        (MieGrunDebyeHotEos(V,T,V0,T0,pHotEos,hotExtraInputs{:}));

    addedThermPressFun = [];
    [Pmod,KTmod,Cvmod,gammod] = calcPressThermAddEos(V(:),T(:),T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);

    %scatter(Pmod,V,50,T,'o')

    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    priorEos = pEos;
    priorcovEos = diag(Inf*ones(size(pEos)));

    %pInitEos = pEos + [0 -15 -.5 0 -.2 .3 0];
    pInitEos = pEos + [0  0 0 0 -.2 .3 0];
    fixFlag = zeros(size(priorEos));
    fixFlag(1)   = 1;
    fixFlag(1:3)   = 1;
    fixFlag(4)   = 1;
    fixFlag(end) = 1;
    %priorcovEos(end,end) = 0;

    PerrTot = 1e-6*ones(size(V));
    P = Pmod;

    %PInit = coldEosFun(V,pInitEos);
    opt = [];
    opt.NfitIter = 1;
    [pfit pfitcov] = fitHotCompressData(pInitEos,fixFlag,T0,...
        NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun,P,V,T,PerrTot,opt);
    relErr = (pfit-pEos)./pEos;

    verifyTrue(testCase,all(abs(relErr)<TOL));
end

%function testFitHotCompressData_credRegion(testCase)
%    TOL = 1e-3;
%
%    T0 = 300;
%    V0 = 162.373;
%    K0 = 258.4;
%    KP0= 4.10;
%
%    Natom = 4*5;
%    Tdeb0 = 940;
%    gam0  = 1.55;
%    q = 1.1;
%
%    pColdEos = [V0 K0 KP0];
%    pHotEos  = [Tdeb0 gam0 q 1.0];
%
%    VCold = V0*linspace(.7,1,9)';
%    TCold = 300*ones(size(VCold));
%    THot = [...
%        [1000:1500:5500]'; ...
%        [1000:1500:5500]'; ...
%        ];
%    VHot  = [...
%        VCold(1)*ones(4,1); ...
%        VCold(9)*ones(4,1); ...
%        ];
%    V = [VCold;VHot];
%    T = [TCold;THot];
%
%
%    coldEosFun = @VinetEos;
%    debyeDerivsFun = @debyePowerLaw;
%    hotExtraInputs = {Natom, debyeDerivsFun};
%    hotEosFun  = @(V,T,V0,T0,pHotEos,hotExtraInputs)...
%        (MieGrunDebyeHotEos(V,T,V0,T0,pHotEos,hotExtraInputs{:}));
%
%    addedThermPressFun = [];
%    [Pmod,KTmod,Cvmod,gammod] = calcPressThermAddEos(V(:),T(:),T0,pColdEos,pHotEos,...
%        coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);
%
%    %scatter(Pmod,V,50,T,'o')
%
%    pEos = [pColdEos pHotEos];
%    NpCold = length(pColdEos);
%
%    priorEos = pEos;
%    priorcovEos = diag(Inf*ones(size(pEos)));
%
%    pInitEos = pEos + [10 -15 -.5 0 -.2 .3 0];
%    fixFlag = zeros(size(priorEos));
%    fixFlag(4)   = 1;
%    fixFlag(end) = 1;
%    %priorcovEos(end,end) = 0;
%
%    PerrTot = 1e-6*ones(size(V));
%    P = Pmod;
%
%    %PInit = coldEosFun(V,pInitEos);
%    opt = [];
%    opt.NfitIter = 1;
%    [pfit pfitcov] = fitHotCompressData(pInitEos,fixFlag,T0,...
%        NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
%        addedThermPressFun,P,V,T,PerrTot,opt);
%    relErr = (pfit-pEos)./pEos;
%
%    verifyTrue(testCase,all(abs(relErr)<TOL));
%end
%    plvl = [0:.01:.99];
%    plot([0;costFunHists],[0:Niter]'/Niter,'k-',invchi2(plvl,dof),plvl,'r-')
%    plot([0;chi2MinHists],[0:Niter]'/Niter,'k-' ,invchi2(plvl,length(V)-dof),plvl,'r-')
%    plot([0;chi2TrueHists],[0:Niter]'/Niter,'k-',invchi2(plvl,length(V)),plvl,'r-')
%
%    indMed = round(.5*Niter);
%
%    chi2draw = chi2rnd(dof,Niter,Ntry);
%    chi2draws = sort(chi2draw);
%    hist(chi2draws(indMed,:),50)
%    costFunHists(indMed)
%
%
%    chi2draw = chi2rnd(length(V)-length(pEos),Niter,Ntry);
%    chi2draws = sort(chi2draw);
%    hist(chi2draws(500,:),50)
%    chi2MinHists(indMed)
%
%    %NOTE: not sure what to do since costFun value at true value is distributed
%    %like a chi square distribution with dof = num parameters - some
%
%    Ndraw = 400;
%    costdraw = zeros(Niter,Ndraw);
%    for(i=1:Niter)
%        ipcov = squeeze(pcovHist(i,:,:));
%        iphess = inv(ipcov);
%        ipfit = pEosHist(i,:);
%        idelpdraw = mvnrnd(ipfit-pEos,ipcov,Ndraw);
%        for(j=1:Ndraw)
%            costdraw(i,j) = idelpdraw(j,:)*iphess*idelpdraw(j,:)';
%        end
%    end
%    costdraws = sort(costdraw);
%
%    plot([0;costdraws(:,1)],[0:Niter]'/Niter,'k-',invchi2(plvl,4),plvl,'r-')
%
%    k=1;plot([0;costdraws(:,k)],[0:Niter]'/Niter,'k-',[0;costFunHists],[0:Niter]'/Niter,'b-',invchi2(plvl,dof),plvl,'r-')
%        mvnrnd(
%
%    k=1;plot([0;chi2draws(:,k)],[0:Niter]'/Niter,'k-',invchi2(plvl,dof),plvl,'r-')
%    round(.68*Niter)
%
%    for(i=1:Ntry)
%
%    Niter
%
%    relErr = (pfit-pEos)./pEos;
%
%    verifyTrue(testCase,all(abs(relErr)<TOL));
