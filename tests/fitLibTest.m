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
    [pColdEos,pHotEos,T0,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun] = getMgPvEos();

    Ncold = 8;
    NhotCyc = 3;
    NstepPerCyc = 4;

    V0 = pColdEos(1);
    [V,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0);

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

    [pColdEos,pHotEos,T0,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun] = getMgPvEos();

    Ncold = 8;
    NhotCyc = 3;
    NstepPerCyc = 4;

    V0 = pColdEos(1);
    [V,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0);

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

function testFitErrModPVT_synth(testCase)
    rndSeed = 42;
    rng(rndSeed)
    TOL = 1e-3;

    [pColdEosMgPv,pHotEosMgPv,T0MgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv] = getMgPvEos();
    [pColdEosMgO,pHotEosMgO,T0MgO,coldEosFunMgO,hotEosFunMgO,...
        hotExtraInputsMgO,addedThermPressFunMgO] = getMgOEos();

    Ncold = 8;
    NhotCyc = 3;
    NstepPerCyc = 4;

    Ncold = 60;
    NhotCyc = 8;
    NstepPerCyc = 20;
    V0MgPv = pColdEosMgPv(1);

    [Vsamp,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0MgPv);

    [Psamp,KTsamp,Cvsamp,gamsamp] = calcPressThermAddEos(Vsamp(:),T(:),...
        T0MgPv,pColdEosMgPv,pHotEosMgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv);

    Vmark = invertPressEos(Psamp,T,pColdEosMgO,pHotEosMgO,T0MgO,coldEosFunMgO,...
        hotEosFunMgO,hotExtraInputsMgO,addedThermPressFunMgO);


    %[Pcalc] = calcPressThermAddEos(Vmark,T,T0MgO,pColdEosMgO,pHotEosMgO,...
    %    coldEosFunMgO,hotEosFunMgO,hotExtraInputsMgO,addedThermPressFunMgO);


    indCold300 = find(T==300);
    % Add errors to V and T
    VsampErr = 0.005*Vsamp;
    VmarkErr = 0.005*Vmark;
    TErr     = 150*ones(size(T));
    TErr(indCold300) = 0;

    ptrueErrMod = [.2 -.2 .0];

    %VmarkObs = Vmark + VmarkErr.*randn(size(VmarkErr));
    %VsampObs = Vsamp + VsampErr.*randn(size(VsampErr));
    %TObs     = T     + TErr.*randn(size(TErr));
    VmarkObs = Vmark + exp(ptrueErrMod(1))*VmarkErr.*randn(size(VmarkErr));
    VsampObs = Vsamp + exp(ptrueErrMod(2))*VsampErr.*randn(size(VsampErr));
    TObs     = T     + exp(ptrueErrMod(3))*TErr.*randn(size(TErr));

    [PmarkObs,KTmarkObs,CvmarkObs,gammarkObs,thmExpmarkObs] = ...
        calcPressThermAddEos(VmarkObs,TObs,T0MgO,pColdEosMgO,pHotEosMgO,...
        coldEosFunMgO,hotEosFunMgO,hotExtraInputsMgO,addedThermPressFunMgO);

    [PsampObs,KTsampObs,CvsampObs,gamsampObs,thmExpsampObs] = ...
        calcPressThermAddEos(VsampObs,TObs,T0MgPv,pColdEosMgPv,pHotEosMgPv,...
        coldEosFunMgPv,hotEosFunMgPv,hotExtraInputsMgPv,addedThermPressFunMgPv);

    dPmarkdVObs    = -KTmarkObs./VmarkObs;
    dPmarkdTObs    = KTmarkObs.*thmExpmarkObs;
    PmarkderivsObs = [dPmarkdVObs,dPmarkdTObs];

    dPsampdVObs    = -KTsampObs./VsampObs;
    dPsampdTObs    = KTsampObs.*thmExpsampObs;
    PsampderivsObs = [dPsampdVObs,dPsampdTObs];
    
    %scatter(PmarkObs,VmarkObs,50,T,'o','filled')
    %scatter(PsampObs,VsampObs,50,T,'o','filled')

    PresidObs = PmarkObs - PsampObs;

    measGrpID = ones(length(PresidObs),1);
    pinitErrMod = [0 0 0];
    priorErrMod = pinitErrMod;
    priorcovErrMod = diag((.3*[1 1 1]).^2);
    opt = [];
    [pfitErrMod, pfitcovErrMod] = fitErrModPVT(pinitErrMod,...
        priorErrMod,priorcovErrMod,PresidObs,PmarkderivsObs,PsampderivsObs,...
        VmarkErr,VsampErr,TErr,measGrpID,opt);
    pfitErrModCredWid = sqrt(diag(pfitcovErrMod)');

    relDev = pfitErrMod./pfitErrModCredWid;

    verifyTrue(testCase,all(abs(relDev)<1),...
        ['Deviation from zero must be less than '...
        'credible region width for errorMod params']);
end



function [pColdEos,pHotEos,T0,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun] = getMgPvEos()
    T0 = 300;

    V0 = 162.373;
    K0 = 258.4;
    KP0= 4.10;

    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyePowerLaw;
    Natom = 4*5;
    Tdeb0 = 940;
    gam0  = 1.55;
    q = 1.1;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 q 1.0];
    
    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    % DOES NOT INCLUDE CORRELATION
    pEosCov = diag([0 1.7 0.07 140 0.09 0.3 0]);

    coldEosFun = @VinetEos;
    hotEosFun = @MieGrunDebyeHotEos;
    debyeDerivsFun = @debyePowerLaw;

    hotExtraInputs = {Natom, debyeDerivsFun};
    addedThermPressFun = [];
end

function [pColdEos,pHotEos,T0,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun] = getMgOEos()

    T0 = 300;
    V0 = 74.698;
    K0 = 160.63;
    KP0= 4.367;

    Natom = 4*2;
    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 a b 1.0];

    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    % DOES NOT INCLUDE CORRELATION
    pEosCov = diag([0 0.18 0.013 13 0.015 0.019 1.1 0]);

    coldEosFun = @VinetEos;
    hotEosFun = @MieGrunDebyeHotEos;
    debyeDerivsFun = @debyeTange;
    
    hotExtraInputs = {Natom, debyeDerivsFun};
    addedThermPressFun = [];

end

function [V,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0)

    ThotBnd = [1000,5000];

    VCold = V0*linspace(.7,1,Ncold)';
    TCold = 300*ones(Ncold,1);
    THotCycle = linspace(ThotBnd(1),ThotBnd(2),NstepPerCyc)';
    THot = repmat(THotCycle(:),NhotCyc,1);

    indHotAnchor = round(linspace(1,Ncold,NhotCyc)');

    VHot = reshape(repmat(VCold(indHotAnchor),1,NstepPerCyc),[],1);

    V = [VCold;VHot];
    T = [TCold;THot];
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
%    hotEosFun  = @MieGrunDebyeHotEos;
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
