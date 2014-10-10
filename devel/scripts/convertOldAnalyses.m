load('../data/working.mat')
datT12M=[pvdatAllR(1).runID pvdatAllR(1).P pvdatAllR(1).Perr pvdatAllR(1).V pvdatAllR(1).Verr pvdatAllR(1).T pvdatAllR(1).Terr  pvdatAllR(1).PerrTot]
datK09M=[pvdatAllR(2).runID pvdatAllR(2).P pvdatAllR(2).Perr pvdatAllR(2).V pvdatAllR(2).Verr pvdatAllR(2).T pvdatAllR(2).Terr  pvdatAllR(2).PerrTot]

sprintf('%1d  %7.3f %4.3f  %6.2f %4.3f  %4.0f %2.0f   %.3f\n',datT12M')
sprintf('%1d  %7.3f %4.3f  %6.2f %4.3f  %4.0f %2.0f   %.3f\n',datK09M')

%Save data in text files via copy paste
% reload data from text files

clear datT12;
datT12.measGrpID = datT12M(:,1);
datT12.P         = datT12M(:,2);
datT12.Perr      = datT12M(:,3);
datT12.V         = datT12M(:,4);
datT12.Verr      = datT12M(:,5);
datT12.T         = datT12M(:,6);
datT12.Terr      = datT12M(:,7);
datT12.PerrTot   = datT12M(:,8);


clear datK09;
datK09.measGrpID = datK09M(:,1);
datK09.P         = datK09M(:,2);
datK09.Perr      = datK09M(:,3);
datK09.V         = datK09M(:,4);
datK09.Verr      = datK09M(:,5);
datK09.T         = datK09M(:,6);
datK09.Terr      = datK09M(:,7);
datK09.PerrTot   = datK09M(:,8);

scatter(datK09.P,datK09.V,50,datK09.T,'o')
scatter(datT12.P,datT12.V,50,datT12.T,'o')

% Fit these data
% Initial fit params taken from Tange 2012
V0 = 162.373;
K0 = 258.4;
KP0= 4.10;

Natom = 4*5;
Tdeb0 = 940;
gam0  = 1.55;
q = 1.1;
Cvfac = 1;

pColdEos = [V0 K0 KP0];
pHotEos  = [Tdeb0 gam0 q Cvfac];

coldEosFun = @VinetEos;
debyeDerivsFun = @debyePowerLaw;
hotExtraInputs = {Natom, debyeDerivsFun};
hotEosFun  = @(V,T,V0,T0,pHotEos,hotExtraInputs)...
    (MieGrunDebyeHotEos(V,T,V0,T0,pHotEos,hotExtraInputs{:}));

addedThermPressFun = [];
%[Pmod,KTmod,Cvmod,gammod] = calcPressThermAddEos(V(:),T(:),T0,pColdEos,pHotEos,...
%    coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);

%scatter(Pmod,V,50,T,'o')

pEos = [pColdEos pHotEos];
NpCold = length(pColdEos);

%priorEos = pEos + [0 5 -.1 100 .1 -.1 0];
priorEos = pEos;
priorcovEos = diag(Inf*ones(size(pEos)));
%priorcovEos(5,5) = 1;
%priorcovEos(6,6) = 1;

pInitEos = priorEos;
fixFlag = zeros(size(priorEos));
fixFlag(1)   = 1;
fixFlag(1:3)   = 1;
%fixFlag(4)   = 1;
fixFlag(end) = 1;


%PInit = coldEosFun(V,pInitEos);
opt = [];
opt.NfitIter = 1;
%opt.robustFit = 1;
[pfitT12 pfitcovT12] = fitHotCompressData(pInitEos,fixFlag,T0,...
    NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
    addedThermPressFun,datT12.P,datT12.V,datT12.T,datT12.PerrTot,opt);
%[pfitT12 pfitcovT12] = fitHotCompressData(pInitEos,fixFlag,T0,...
%    NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
%    addedThermPressFun,datT12.P,datT12.V,datT12.T,datT12.Perr,opt);

%pfitT12
%sqrt(diag(pfitcovT12))

pfitT12(1:3)
pInitEos(1:3)
pfitT12(4)
pInitEos(4)
pfitT12(5:6)
pInitEos(5:6)

