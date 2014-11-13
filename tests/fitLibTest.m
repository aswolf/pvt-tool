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
    [V,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0,[],[]);

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
    [V,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0,[],[]);

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
    TOLWID = 0.3;

    [pColdEosMgPv,pHotEosMgPv,T0MgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv] = getMgPvEos();
    [pColdEosMgO,pHotEosMgO,T0MgO,coldEosFunMgO,hotEosFunMgO,...
        hotExtraInputsMgO,addedThermPressFunMgO] = getMgOEos();

    Ndraw = 30;

    Ncold = 8;
    NhotCyc = 3;
    NstepPerCyc = 4;

    Ncold = 60;
    NhotCyc = 8;
    NstepPerCyc = 20;
    V0MgPv = pColdEosMgPv(1);

    [Vsamp,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0MgPv,[],[]);

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
    relDev = zeros(Ndraw,2);
    for(i=1:Ndraw)
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

        measGrpID = 1;
        pinitErrMod = [0 0];
        priorErrMod = pinitErrMod;
        priorcovErrMod = diag((.3*[1 1]).^2);
        opt = [];
        [pfitErrMod, pfitcovErrMod] = fitErrModPVT(measGrpID,pinitErrMod,...
            priorErrMod,priorcovErrMod,PresidObs,PmarkderivsObs,PsampderivsObs,...
            VmarkErr,VsampErr,TErr,opt);
        pfitErrModCredWid = sqrt(diag(pfitcovErrMod)');

        irelDev = pfitErrMod./pfitErrModCredWid;
        relDev(i,:) = irelDev;
    end

    verifyTrue(testCase,abs(mean(relDev(:,1)))/std(relDev(:,1)) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    verifyTrue(testCase,abs(log(std(relDev(:,1)))) < TOLWID ,...
        ['Distribution of normalized errMod param for avg magnitude must '...
        'have width close to truth (ie. within TOLWID of 1)']);
end

function testFitErrModPVT_synthMultiMeasGrp(testCase)
    rndSeed = 42;
    rng(rndSeed)
    TOLWID = 0.3;

    [pColdEosMgPv,pHotEosMgPv,T0MgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv] = getMgPvEos();
    [pColdEosMgO,pHotEosMgO,T0MgO,coldEosFunMgO,hotEosFunMgO,...
        hotExtraInputsMgO,addedThermPressFunMgO] = getMgOEos();

    Ndraw = 30;

    Ncold = 8;
    NhotCyc = 3;
    NstepPerCyc = 4;

    Ncold = 60;
    NhotCyc = 8;
    NstepPerCyc = 20;
    V0MgPv = pColdEosMgPv(1);

    [Vsamp,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0MgPv,[],[]);

    [Psamp,KTsamp,Cvsamp,gamsamp] = calcPressThermAddEos(Vsamp(:),T(:),...
        T0MgPv,pColdEosMgPv,pHotEosMgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv);

    Vmark = invertPressEos(Psamp,T,pColdEosMgO,pHotEosMgO,T0MgO,coldEosFunMgO,...
        hotEosFunMgO,hotExtraInputsMgO,addedThermPressFunMgO);


    %[Pcalc] = calcPressThermAddEos(Vmark,T,T0MgO,pColdEosMgO,pHotEosMgO,...
    %    coldEosFunMgO,hotEosFunMgO,hotExtraInputsMgO,addedThermPressFunMgO);

    indCold300 = find(T==300);
    % Add errors to V and T
    VsampErr = 0.003*Vsamp;
    VmarkErr = 0.003*Vmark;
    TErr     = 80*ones(size(T));
    % Zero out errors for ambient temp meas
    TErr(indCold300) = 0;

    Ndat = length(Vsamp);

    % keyboard;
    % Randomly set half of data to second group
    measGrpInd = ones(Ndat,1);
    %measGrpInd(rand(Ndat,1)<0.5) = 2;
    measGrpInd(T>310) = 2;
    measGrpID = cellstr(num2str(measGrpInd));

    % ptrueErrMod = [.2 -.2 .0];
    %               grp1      grp2
    %              dV  dT    dV  dT 
    ptrueErrMod = [.2 -.2; -0.1 -.25];
    ptrueErrMod = [.2 0; -0.1 .4];

    %VmarkObs = Vmark + VmarkErr.*randn(size(VmarkErr));
    %VsampObs = Vsamp + VsampErr.*randn(size(VsampErr));
    %TObs     = T     + TErr.*randn(size(TErr));
    relDev  = zeros(Ndraw,2,2);
    fitdraw = zeros(Ndraw,2,2);
    NGrp1 = sum(measGrpInd==1);
    NGrp2 = sum(measGrpInd==2);
    VmarkDev = zeros(Ndat,1);
    VsampDev = zeros(Ndat,1);
    TDev     = zeros(Ndat,1);
    for(i=1:Ndraw)
        VmarkDev(measGrpInd==1)=exp(ptrueErrMod(1,1))*VmarkErr(measGrpInd==1).*randn(NGrp1,1);
        VmarkDev(measGrpInd==2)=exp(ptrueErrMod(2,1))*VmarkErr(measGrpInd==2).*randn(NGrp2,1);
        VsampDev(measGrpInd==1)=exp(ptrueErrMod(1,1))*VsampErr(measGrpInd==1).*randn(NGrp1,1);
        VsampDev(measGrpInd==2)=exp(ptrueErrMod(2,1))*VsampErr(measGrpInd==2).*randn(NGrp2,1);
        TDev(measGrpInd==1)=exp(ptrueErrMod(1,2))*TErr(measGrpInd==1).*randn(NGrp1,1);
        TDev(measGrpInd==2)=exp(ptrueErrMod(2,2))*TErr(measGrpInd==2).*randn(NGrp2,1);

        VmarkObs = Vmark + VmarkDev;
        VsampObs = Vsamp + VsampDev;
        TObs     = T     + TDev;
        %log(std(TDev(measGrpInd==2)')/TErr(end))       

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

        PresidDerivs = [PmarkderivsObs(:,1) PsampderivsObs(:,1) ...
            (PmarkderivsObs(:,2)-PsampderivsObs(:,2))];

        pinitErrModList = [0 0; 0 0];
        priorErrModList = pinitErrModList;
        priorcovErrModList = zeros(2,2,2);
        priorcovErrModList(1,:,:) = diag(([0.3 0.3]).^2);
        priorcovErrModList(2,:,:) = diag(([0.3 0.3]).^2);
        opt = [];
        [pfitErrModList, pfitcovErrModList, nLogPFunList] = fitErrModPVT(measGrpID,pinitErrModList,...
            priorErrModList,priorcovErrModList,PresidObs,PmarkderivsObs,PsampderivsObs,...
            VmarkErr,VsampErr,TErr,opt);
        pfitErrModCredWidList(1,:) = sqrt(diag(squeeze(pfitcovErrModList(1,:,:)))');
        pfitErrModCredWidList(2,:) = sqrt(diag(squeeze(pfitcovErrModList(2,:,:)))');

        fitdraw(i,:,:) = pfitErrModList;
        irelDev = (pfitErrModList-ptrueErrMod(:,1:2)) ./pfitErrModCredWidList;
        relDev(i,:,:) = irelDev;
    end
    %ptrueErrMod
    %hist(fitdraw(:,1,1),8)
    %hist(fitdraw(:,2,1),8)
    %hist(fitdraw(:,1,2),8)
    %hist(fitdraw(:,2,2),8)

    %hist(relDev(:,:,1),8)
    %hist(relDev(:,:,2),8)
    %hist(relDev(:,1,2),8)

    verifyTrue(testCase,false ,...
        ['MultiMeasGrp tested by hand but not automated yet']);
    verifyTrue(testCase,abs(mean(relDev(:,1)))/std(relDev(:,1)) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    verifyTrue(testCase,abs(log(std(relDev(:,1)))) < TOLWID ,...
        ['Distribution of normalized errMod param for avg magnitude must '...
        'have width close to truth (ie. within TOLWID of 1)']);
end

function testFitErrModPVT_synthMultiMeasGrpSoftPmark(testCase)
    rndSeed = 42;
    rng(rndSeed)
    TOLWID = 0.3;

    [pColdEosMgPv,pHotEosMgPv,T0MgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv] = getMgPvEos();
    [pColdEosNe,pHotEosNe,T0Ne,coldEosFunNe,hotEosFunNe,...
        hotExtraInputsNe,addedThermPressFunNe] = getNeEos();

    Ndraw = 30;

    Ncold = 8;
    NhotCyc = 3;
    NstepPerCyc = 4;

    Ncold = 60;
    NhotCyc = 8;
    NstepPerCyc = 20;
    V0MgPv = pColdEosMgPv(1);

    [Vsamp,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0MgPv,...
        [.75 .9],[1000 3000]);

    [Psamp,KTsamp,Cvsamp,gamsamp] = calcPressThermAddEos(Vsamp(:),T(:),...
        T0MgPv,pColdEosMgPv,pHotEosMgPv,coldEosFunMgPv,hotEosFunMgPv,...
        hotExtraInputsMgPv,addedThermPressFunMgPv);

    Vmark = invertPressEos(Psamp,T,pColdEosNe,pHotEosNe,T0Ne,coldEosFunNe,...
        hotEosFunNe,hotExtraInputsNe,addedThermPressFunNe);

    [Pcalc] = calcPressThermAddEos(Vmark,T,T0Ne,pColdEosNe,pHotEosNe,...
        coldEosFunNe,hotEosFunNe,hotExtraInputsNe,addedThermPressFunNe);
    
    indCold300 = find(T==300);
    % Add errors to V and T
    VsampErr = 0.003*Vsamp;
    VmarkErr = 0.003*Vmark;
    TErr     = 80*ones(size(T));
    % Zero out errors for ambient temp meas
    TErr(indCold300) = 0;

    Ndat = length(Vsamp);

    % keyboard;
    % Randomly set half of data to second group
    measGrpInd = ones(Ndat,1);
    %measGrpInd(rand(Ndat,1)<0.5) = 2;
    measGrpInd(T>310) = 2;
    measGrpID = cellstr(num2str(measGrpInd));

    % ptrueErrMod = [.2 -.2 .0];
    %               grp1      grp2
    %              dV  dT    dV  dT 
    ptrueErrMod = [.2 -.2; -0.1 -.25];
    ptrueErrMod = [.2 0; -0.1 -.2];

    %VmarkObs = Vmark + VmarkErr.*randn(size(VmarkErr));
    %VsampObs = Vsamp + VsampErr.*randn(size(VsampErr));
    %TObs     = T     + TErr.*randn(size(TErr));
    relDev  = zeros(Ndraw,2,2);
    fitdraw = zeros(Ndraw,2,2);
    PErrTermMagList = zeros(Ndraw,2,3);
    NGrp1 = sum(measGrpInd==1);
    NGrp2 = sum(measGrpInd==2);
    VmarkDev = zeros(Ndat,1);
    VsampDev = zeros(Ndat,1);
    TDev     = zeros(Ndat,1);
    for(i=1:Ndraw)
        VmarkDev(measGrpInd==1)=exp(ptrueErrMod(1,1))*VmarkErr(measGrpInd==1).*randn(NGrp1,1);
        VmarkDev(measGrpInd==2)=exp(ptrueErrMod(2,1))*VmarkErr(measGrpInd==2).*randn(NGrp2,1);
        VsampDev(measGrpInd==1)=exp(ptrueErrMod(1,1))*VsampErr(measGrpInd==1).*randn(NGrp1,1);
        VsampDev(measGrpInd==2)=exp(ptrueErrMod(2,1))*VsampErr(measGrpInd==2).*randn(NGrp2,1);
        TDev(measGrpInd==1)=exp(ptrueErrMod(1,2))*TErr(measGrpInd==1).*randn(NGrp1,1);
        TDev(measGrpInd==2)=exp(ptrueErrMod(2,2))*TErr(measGrpInd==2).*randn(NGrp2,1);

        VmarkObs = Vmark + VmarkDev;
        VsampObs = Vsamp + VsampDev;
        TObs     = T     + TDev;
        %log(std(TDev(measGrpInd==2)')/TErr(end))       

        [PmarkObs,KTmarkObs,CvmarkObs,gammarkObs,thmExpmarkObs] = ...
            calcPressThermAddEos(VmarkObs,TObs,T0Ne,pColdEosNe,pHotEosNe,...
            coldEosFunNe,hotEosFunNe,hotExtraInputsNe,addedThermPressFunNe);

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

        PresidDerivs = [PmarkderivsObs(:,1) PsampderivsObs(:,1) ...
            (PmarkderivsObs(:,2)-PsampderivsObs(:,2))];
        iPErrTermMag = abs(PresidDerivs).*[VmarkErr VsampErr TErr];
        PErrTermMagList(i,1,:)  = median(iPErrTermMag(measGrpInd==1,:));
        PErrTermMagList(i,2,:)  = median(iPErrTermMag(measGrpInd==2,:));

        pinitErrModList = [0 0; 0 0];
        priorErrModList = pinitErrModList;
        priorcovErrModList = zeros(2,2,2);
        priorcovErrModList(1,:,:) = diag(([0.3 0.3]).^2);
        priorcovErrModList(2,:,:) = diag(([0.3 0.3]).^2);
        opt = [];
        [pfitErrModList, pfitcovErrModList, nLogPFunList] = fitErrModPVT(measGrpID,pinitErrModList,...
            priorErrModList,priorcovErrModList,PresidObs,PmarkderivsObs,PsampderivsObs,...
            VmarkErr,VsampErr,TErr,opt);
        pfitErrModCredWidList(1,:) = sqrt(diag(squeeze(pfitcovErrModList(1,:,:)))');
        pfitErrModCredWidList(2,:) = sqrt(diag(squeeze(pfitcovErrModList(2,:,:)))');

        fitdraw(i,:,:) = pfitErrModList;
        irelDev = (pfitErrModList-ptrueErrMod(:,1:2)) ./pfitErrModCredWidList;
        relDev(i,:,:) = irelDev;
    end
    PErrTermMag(1,:) = median(squeeze(PErrTermMagList(:,1,:)));
    PErrTermMag(2,:) = median(squeeze(PErrTermMagList(:,2,:)));
    %PErrTermMag
    %ptrueErrMod
    %hist(fitdraw(:,1,1),8)
    %hist(fitdraw(:,2,1),8)
    %hist(fitdraw(:,1,2),8)
    %hist(fitdraw(:,2,2),8)

    %hist(relDev(:,:,1),8)
    %hist(relDev(:,:,2),8)
    %hist(relDev(:,2,2),8)

    verifyTrue(testCase,false ,...
        ['MultiMeasGrp Soft Pmark tested by hand but not automated yet']);
    verifyTrue(testCase,abs(mean(relDev(:,1)))/std(relDev(:,1)) < 1 ,...
        ['Distribution of normalized errMod params must have small systematic '...
        'deviation from truth (ie. within 1 sigma of zero)']);
    verifyTrue(testCase,abs(log(std(relDev(:,1)))) < TOLWID ,...
        ['Distribution of normalized errMod param for avg magnitude must '...
        'have width close to truth (ie. within TOLWID of 1)']);
end

function [pColdEos,pHotEos,T0,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun] = getNeEos()

    T0 = 300;
    V0 = 22.234;
    K0 = 1.070;
    KP0= 8.40;

    Natom = 4;
    Tdeb0 = 75.1;
    gam0  = 2.442;
    q = 0.97;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 q 1.0];
    
    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    % DOES NOT INCLUDE CORRELATION
    pEosCov = diag([0 0.016 0.03 0 0 0 0]);

    coldEosFun = @VinetEos;
    hotEosFun = @MieGrunDebyeHotEos;
    debyeDerivsFun = @debyePowerLaw;

    hotExtraInputs = {Natom, debyeDerivsFun};
    addedThermPressFun = [];

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

function [V,T] = makeIsobarPVTdata(Ncold,NhotCyc,NstepPerCyc,V0,VfacBnd,ThotBnd)
    if(isempty(VfacBnd))
        VfacBnd = [.7,1];
    end
    if(isempty(ThotBnd))
        ThotBnd = [1000, 5000];
    end

    VCold = V0*linspace(VfacBnd(1), VfacBnd(2),Ncold)';
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
