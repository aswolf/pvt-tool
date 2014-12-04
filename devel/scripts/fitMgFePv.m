
eosMgPvT12 = getEos_MgPvTange2012();
eosMgOMod = getEos_MgOTange2009();
eosNeMod = getEos_NeDewaele2008();

dataDir = '~/Documents/code/MATLAB/pvt-tool/devel/data';

loadScriptNmK09 = fullfile(dataDir,'loadPVT_Katsura2009.m');
markEos = eosMgOMod;
MgPvPVTdataK09 = loadPVTdata(loadScriptNmK09,markEos);

loadScriptNmT12 = fullfile(dataDir,'loadPVT_Tange2012.m');
markEos = eosMgOMod;
MgPvPVTdataT12 = loadPVTdata(loadScriptNmT12,markEos);

loadScriptNmW14 = fullfile(dataDir,'loadPVT_Wolf2014.m');
markEos = eosNeMod;
MgFePvPVTdataW14 = loadPVTdata(loadScriptNmW14,markEos);



filename = fullfile(dataDir,'fit_Wolf2014_lsq.in');
keywords = {'fileName','name','measGrpIDLbl','errMode','markLbl'};
file_struct=parse_file(filename,keywords)

fid = fopen(loadScriptNmW14,'r')
C = textscan(fid, '%s','Delimiter','');
fclose(fid)
C = C{:};
C = cellfun(@strtrim,C,'UniformOutput',false);


% NOTE we need Pt EOS
%loadScriptNmF98 = fullfile(dataDir,'loadPVT_Fiquet1998.m');
%markEos = eosPtMod;
%MgPvPVTdataF98 = loadPVTdata(loadScriptNmF98,markEos);

%loadScriptNmL08 = fullfile(dataDir,'loadPVT_Lundin2008_00.m');
%markEos = [];
%MgPvPVTdataL08 = loadPVTdata(loadScriptNmL08,markEos);

%MgPvPVTdata = getPVTdata_MgPvTange2012('mark');
% FePvPVTdata = getPVTdata_MgFePvWolf2014();
% MgPvPVTdataF98 = getPVTdata_MgPvFiquet1998();
% load Tange Eos for MgPv and MgO for initial guess and marker EOS

% Perform multiple tests with different priors
% Fit Wolf Data: Prior all Tange values
name1W14 = 'MgFePv: Refit all params with good priors (Wolf2014 data)';
fixFlag1W14 = [0 0 0 1 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [163.2];
pEosPrior(4) = 1000;
pEosPrior(5:6) = [1 1];
eosMgFePvPrior1W14 = eosMgPvT12;
eosMgFePvPrior1W14.pEos = pEosPrior;
eosMgFePvPrior1W14.pEosCov = zeros(7);
eosMgFePvPrior1W14.pEosCov(1,1) = 0.2^2;
eosMgFePvPrior1W14.pEosCov(2,2) = Inf;
eosMgFePvPrior1W14.pEosCov(3,3) = Inf;
eosMgFePvPrior1W14.pEosCov(4,4) = 0^2;
eosMgFePvPrior1W14.pEosCov(5,5) =   1^2;
eosMgFePvPrior1W14.pEosCov(6,6) =   1^2;

% Initialize and perform fit using least squares approach
opt = [];
PVTeval1W14 = initPVTeval(name1W14,MgFePvPVTdataW14,eosMgFePvPrior1W14,opt);
PVTeval1W14 = fitSampModPVTeval(PVTeval1W14,[],fixFlag1W14);
% Review residuals to check if error model is needed
hist((PVTeval1W14.Psamp-PVTeval1W14.PVTdataList.Pmark)./PVTeval1W14.PVTdataList.PErrTot,15)
std((PVTeval1W14.Psamp-PVTeval1W14.PVTdataList.Pmark)./PVTeval1W14.PVTdataList.PErrTot)
PVTeval1W14 = fitErrModPVTeval(PVTeval1W14,[]);
hist((PVTeval1W14.Psamp-PVTeval1W14.PVTdataList.Pmark)./PVTeval1W14.PVTdataList.PErrTot,15)
std((PVTeval1W14.Psamp-PVTeval1W14.PVTdataList.Pmark)./PVTeval1W14.PVTdataList.PErrTot)
% Refit with updated error model
PVTeval1W14 = fitSampModPVTeval(PVTeval1W14,[],fixFlag1W14);
fileNm = '../data/Wolf2014ErrTotFit1.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval1W14.PVTdataList,fileNm);
fileNm = '../data/Wolf2014Fit1.md';
eosModList = [PVTeval1W14.sampEosFit,PVTeval1W14.sampEosPrior];
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'decimal');
%[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'paren');


% Plot results
Presid = PVTeval1W14.Psamp - PVTeval1W14.PVTdataList.Pmark;
[Pred] = evalPressEos([],PVTeval1W14.sampEosFit,MgFePvPVTdataW14.V,300);
Vbnds = [min(MgFePvPVTdataW14.V),max(MgFePvPVTdataW14.V)];
Vrng = diff(Vbnds);
%Vbnds = [Vbnds(1)-.1*Vrng, Vbnds(2)+.1*Vrng];
Vbnds = [Vbnds(1)-.1*Vrng, PVTeval1W14.sampEosFit.pEos(1)];
Vmod = linspace(Vbnds(1),Vbnds(2),100)';
[PredMod] = evalPressEos([],PVTeval1W14.sampEosFit,Vmod,300);
indHot = find(MgFePvPVTdataW14.T>310);

% Plot Reduced Isotherm
clf;
plot(PredMod,Vmod,'b-');
hold on;
caxis([300 2500]);
hh = ploterrcol(Pred-Presid, MgFePvPVTdataW14.V,...
    PVTeval1W14.PVTdataList.PErrTot, [],MgFePvPVTdataW14.T, 'o','hhx',0);
set(hh,'LineWidth',2,'MarkerSize',5)
hold off;
xlim([30 130])
colorbar;
xlabel('Pressure [GPa]');
ylabel('Volume [A^3]');
title('13%Fe MgPv - Reduced 300 K Isotherm')


% Plot hist
x = linspace(-3,3,100)';
[nhist, xhist] = hist(Presid./PVTeval1W14.PVTdataList.PErrTot,20);
hist(Presid./PVTeval1W14.PVTdataList.PErrTot,20);
hold on;
plot(x, length(Presid)*(xhist(2)-xhist(1))/sqrt(2*pi)*exp(-0.5*x.^2),'r-')
hold off;
xlabel('Relative Residual')
title('13%Fe MgPv')

%scatter(MgFePvPVTdataW14.Pmark, MgFePvPVTdataW14.V, 50, MgFePvPVTdataW14.T, 'o', 'filled')

%clf;
%plot(PredMod,Vmod,'b-');
%hold on;
%scatter(MgFePvPVTdataW14.Pmark, MgFePvPVTdataW14.V, 50, MgFePvPVTdataW14.T, 'x')
%scatter(flipud(Pred-Presid), flipud(MgFePvPVTdataW14.V), 50, flipud(MgFePvPVTdataW14.T), 'o','filled')
%hold off;
%xlim([0 140])




% Perform multiple tests with different priors
% Fit Wolf Data: Prior all Tange values
name1W14R = 'MgFePv: Refit all params with good priors (Wolf2014 data)';
fixFlag1W14R = [0 0 0 1 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [163.2];
pEosPrior(4) = 1000;
pEosPrior(5:6) = [1 1];
eosMgFePvPrior1W14R = eosMgPvT12;
eosMgFePvPrior1W14R.pEos = pEosPrior;
eosMgFePvPrior1W14R.pEosCov = zeros(7);
eosMgFePvPrior1W14R.pEosCov(1,1) = 0.2^2;
eosMgFePvPrior1W14R.pEosCov(2,2) = Inf;
eosMgFePvPrior1W14R.pEosCov(3,3) = Inf;
eosMgFePvPrior1W14R.pEosCov(4,4) = 0^2;
eosMgFePvPrior1W14R.pEosCov(5,5) =   1^2;
eosMgFePvPrior1W14R.pEosCov(6,6) =   1^2;


% Initialize and perform fit
optR.robustFit = true;
% Initialize and perform fit using least squares approach
PVTeval1W14R = initPVTeval(name1W14R,MgFePvPVTdataW14,eosMgFePvPrior1W14R,optR);
PVTeval1W14R = fitSampModPVTeval(PVTeval1W14R,[],fixFlag1W14R);
% Review residuals to check if error model is needed
hist((PVTeval1W14R.Psamp-PVTeval1W14R.PVTdataList.Pmark)./PVTeval1W14R.PVTdataList.PErrTot,15)
std((PVTeval1W14R.Psamp-PVTeval1W14R.PVTdataList.Pmark)./PVTeval1W14R.PVTdataList.PErrTot)
PVTeval1W14R = fitErrModPVTeval(PVTeval1W14R,[]);
hist((PVTeval1W14R.Psamp-PVTeval1W14R.PVTdataList.Pmark)./PVTeval1W14R.PVTdataList.PErrTot,15)
std((PVTeval1W14R.Psamp-PVTeval1W14R.PVTdataList.Pmark)./PVTeval1W14R.PVTdataList.PErrTot)
% Refit with updated error model
PVTeval1W14R = fitSampModPVTeval(PVTeval1W14R,[],fixFlag1W14R);
fileNm = '../data/Wolf2014ErrTotFit1R.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval1W14R.PVTdataList,fileNm);
fileNm = '../data/Wolf2014Fit1R.md';
eosModList = [PVTeval1W14R.sampEosFit,PVTeval1W14R.sampEosPrior];
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'decimal');
%[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'paren');

% Perform multiple tests with different priors
% Fit Wolf Data: Prior all Tange values
name2W14R = 'MgFePv: Refit all params with good priors (Wolf2014 data)';
fixFlag2W14R = [0 0 0 0 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [163.2];
pEosPrior(4) = 1000;
pEosPrior(5:6) = [1 1];
eosMgFePvPrior2W14R = eosMgPvT12;
eosMgFePvPrior2W14R.pEos = pEosPrior;
eosMgFePvPrior2W14R.pEosCov = zeros(7);
eosMgFePvPrior2W14R.pEosCov(1,1) = 0.2^2;
eosMgFePvPrior2W14R.pEosCov(2,2) = Inf;
eosMgFePvPrior2W14R.pEosCov(3,3) = Inf;
eosMgFePvPrior2W14R.pEosCov(4,4) = 100^2;
eosMgFePvPrior2W14R.pEosCov(5,5) =   1^2;
eosMgFePvPrior2W14R.pEosCov(6,6) =   1^2;


% Initialize and perform fit
optR.robustFit = true;
% Initialize and perform fit using least squares approach
PVTeval2W14R = initPVTeval(name2W14R,MgFePvPVTdataW14,eosMgFePvPrior2W14R,optR);
PVTeval2W14R = fitSampModPVTeval(PVTeval2W14R,[],fixFlag2W14R);
% Review residuals to check if error model is needed
hist((PVTeval2W14R.Psamp-PVTeval2W14R.PVTdataList.Pmark)./PVTeval2W14R.PVTdataList.PErrTot,15)
std((PVTeval2W14R.Psamp-PVTeval2W14R.PVTdataList.Pmark)./PVTeval2W14R.PVTdataList.PErrTot)
PVTeval2W14R = fitErrModPVTeval(PVTeval2W14R,[]);
hist((PVTeval2W14R.Psamp-PVTeval2W14R.PVTdataList.Pmark)./PVTeval2W14R.PVTdataList.PErrTot,15)
std((PVTeval2W14R.Psamp-PVTeval2W14R.PVTdataList.Pmark)./PVTeval2W14R.PVTdataList.PErrTot)
% Refit with updated error model
PVTeval2W14R = fitSampModPVTeval(PVTeval2W14R,[],fixFlag2W14R);
fileNm = '../data/Wolf2014ErrTotR.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval2W14R.PVTdataList,fileNm);
fileNm = '../data/Wolf2014FitR.md';
eosModList = [PVTeval2W14R.sampEosFit,PVTeval2W14R.sampEosPrior];
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'decimal');
%[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'paren');



% Perform multiple tests with different priors
% Fit Wolf Data: Prior all Tange values
name2W14 = 'MgFePv: Refit all params with good priors (Wolf2014 data)';
fixFlag2W14 = [0 0 0 1 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [163.2];
pEosPrior(4) = 1100;
pEosPrior(5:6) = [1 1];
eosMgFePvPrior2W14 = eosMgPvT12;
eosMgFePvPrior2W14.pEos = pEosPrior;
eosMgFePvPrior2W14.pEosCov = zeros(7);
eosMgFePvPrior2W14.pEosCov(1,1) = 0.2^2;
eosMgFePvPrior2W14.pEosCov(2,2) = Inf;
eosMgFePvPrior2W14.pEosCov(3,3) = Inf;
eosMgFePvPrior2W14.pEosCov(4,4) = 0^2;
eosMgFePvPrior2W14.pEosCov(5,5) =   1^2;
eosMgFePvPrior2W14.pEosCov(6,6) =   1^2;

% Initialize and perform fit using least squares approach
PVTeval2W14 = initPVTeval(name2W14,MgFePvPVTdataW14,eosMgFePvPrior2W14,opt);
PVTeval2W14 = fitSampModPVTeval(PVTeval2W14,[],fixFlag2W14);
% Review residuals to check if error model is needed
hist((PVTeval2W14.Psamp-PVTeval2W14.PVTdataList.Pmark)./PVTeval2W14.PVTdataList.PErrTot,15)
std((PVTeval2W14.Psamp-PVTeval2W14.PVTdataList.Pmark)./PVTeval2W14.PVTdataList.PErrTot)
PVTeval2W14 = fitErrModPVTeval(PVTeval2W14,[]);
hist((PVTeval2W14.Psamp-PVTeval2W14.PVTdataList.Pmark)./PVTeval2W14.PVTdataList.PErrTot,15)
std((PVTeval2W14.Psamp-PVTeval2W14.PVTdataList.Pmark)./PVTeval2W14.PVTdataList.PErrTot)

% Refit with updated error model
PVTeval2W14 = fitSampModPVTeval(PVTeval2W14,[],fixFlag2W14);
fileNm = '../data/Wolf2014ErrTotFit2.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval2W14.PVTdataList,fileNm);
fileNm = '../data/Wolf2014Fit2.md';
eosModList = [PVTeval2W14.sampEosFit,PVTeval2W14.sampEosPrior];
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'decimal');


% Test1: Prior all Tange values
name1 = 'MgPv: Refit all params with good priors (Tange2012 data)';
fixFlag1 = [0 0 0 0 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [162.5];
pEosPrior(4) = 1100;
pEosPrior(5:6) = [1 1];
eosMgPvPrior1 = eosMgPvT12;
eosMgPvPrior1.pEos = pEosPrior;
eosMgPvPrior1.pEosCov = zeros(7);
eosMgPvPrior1.pEosCov(1,1) = 0.2^2;
eosMgPvPrior1.pEosCov(2,2) = Inf;
eosMgPvPrior1.pEosCov(3,3) = Inf;
eosMgPvPrior1.pEosCov(4,4) = 200^2;
eosMgPvPrior1.pEosCov(5,5) =   1^2;
eosMgPvPrior1.pEosCov(6,6) =   1^2;
% Initialize and perform fit
PVTeval1T12 = initPVTeval(name1,MgPvPVTdataT12,eosMgPvPrior1,opt);
PVTeval1T12 = fitSampModPVTeval(PVTeval1T12,[],fixFlag1);
% Review residuals to check if error model is needed
hist((PVTeval1T12.Psamp-PVTeval1T12.PVTdataList.Pmark)./PVTeval1T12.PVTdataList.PErrTot,15)
std((PVTeval1T12.Psamp-PVTeval1T12.PVTdataList.Pmark)./PVTeval1T12.PVTdataList.PErrTot)
PVTeval1T12 = fitErrModPVTeval(PVTeval1T12,[]);
hist((PVTeval1T12.Psamp-PVTeval1T12.PVTdataList.Pmark)./PVTeval1T12.PVTdataList.PErrTot,15)
std((PVTeval1T12.Psamp-PVTeval1T12.PVTdataList.Pmark)./PVTeval1T12.PVTdataList.PErrTot)
% Refit with updated error model
PVTeval1T12 = fitSampModPVTeval(PVTeval1T12,[],fixFlag1);
fileNm = '../data/Tange2012ErrTot.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval1T12.PVTdataList,fileNm);
fileNm = '../data/Tange2012Fit.md';
eosModList = [PVTeval1T12.sampEosFit,PVTeval1T12.sampEosPrior];
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'decimal');
%[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'paren');



% Plot results
Presid = PVTeval1T12.Psamp - PVTeval1T12.PVTdataList.Pmark;
[Pred] = evalPressEos([],PVTeval1T12.sampEosFit,MgPvPVTdataT12.V,300);
Vbnds = [min(MgPvPVTdataT12.V),max(MgPvPVTdataT12.V)];
Vrng = diff(Vbnds);
%Vbnds = [Vbnds(1)-.1*Vrng, Vbnds(2)+.1*Vrng];
Vbnds = [Vbnds(1)-.1*Vrng, PVTeval1T12.sampEosFit.pEos(1)];
Vmod = linspace(Vbnds(1),Vbnds(2),100)';
[PredMod] = evalPressEos([],PVTeval1T12.sampEosFit,Vmod,300);
indHot = find(MgPvPVTdataT12.T>310);

% Plot Reduced Isotherm
clf;
plot(PredMod,Vmod,'b-');
hold on;
caxis([300 2500]);
hh = ploterrcol(Pred-Presid, MgPvPVTdataT12.V,...
    PVTeval1T12.PVTdataList.PErrTot, [],MgPvPVTdataT12.T, 'o','hhx',0);
set(hh,'LineWidth',2,'MarkerSize',5)
hold off;
xlim([20 100])
colorbar;
xlabel('Pressure [GPa]');
ylabel('Volume [A^3]');
title('MgPv (Tange2012) - Reduced 300 K Isotherm')

% Plot hist
x = linspace(-3,3,100)';
[nhist, xhist] = hist(Presid./PVTeval1T12.PVTdataList.PErrTot,20);
hist(Presid./PVTeval1T12.PVTdataList.PErrTot,20);
hold on;
plot(x, length(Presid)*(xhist(2)-xhist(1))/sqrt(2*pi)*exp(-0.5*x.^2),'r-')
hold off;
xlabel('Relative Residual')
title('MgPv (Tange2012)')


% Initialize and perform fit

%

% Test2: Prior all Tange values
name2 = 'MgPv: Refit params with good priors, fix Tdeb (Tange2012 data)';
fixFlag2 = [0 0 0 1 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [162.5];
pEosPrior(4) = 1100;
pEosPrior(5:6) = [1 1];
eosMgPvPrior2 = eosMgPvT12;
eosMgPvPrior2.pEos = pEosPrior;
eosMgPvPrior2.pEosCov = zeros(7);
eosMgPvPrior2.pEosCov(1,1) = 0.2^2;
eosMgPvPrior2.pEosCov(2,2) = Inf;
eosMgPvPrior2.pEosCov(3,3) = Inf;
eosMgPvPrior2.pEosCov(4,4) =   0;
eosMgPvPrior2.pEosCov(5,5) =   1^2;
eosMgPvPrior2.pEosCov(6,6) =   1^2;

PVTdataFix = PVTdata;
% Test3: Prior all Tange values
name3 = 'MgPv: Refit params with good priors, (Tange2012 data)';
fixFlag3 = [0 0 0 1 0 0 1];
pEosPrior     = eosMgPvT12.pEos;
pEosPrior(1) = [162.5];
pEosPrior(4) = 1000;
pEosPrior(5:6) = [1 1];
eosMgPvPrior3 = eosMgPvT12;
eosMgPvPrior3.pEos = pEosPrior;
eosMgPvPrior3.pEosCov = zeros(7);
eosMgPvPrior3.pEosCov(1,1) = 0.2^2;
eosMgPvPrior3.pEosCov(2,2) = Inf;
eosMgPvPrior3.pEosCov(3,3) = Inf;
eosMgPvPrior3.pEosCov(4,4) =   0;
eosMgPvPrior3.pEosCov(5,5) =   1^2;
eosMgPvPrior3.pEosCov(6,6) =   1^2;

% Initialize and perform fit
PVTeval3 = initPVTeval(name3,PVTeval1.PVTdataList,eosMgPvPrior3,opt);
PVTeval3 = fitSampModPVTeval(PVTeval3,[],fixFlag3);
PVTeval3.sampEosFit
% Review residuals to check if error model is needed
hist((PVTeval1.Psamp-PVTeval1.PVTdataList.Pmark)./PVTeval3.PVTdataList.PErrTot,15)
std((PVTeval1.Psamp-PVTeval1.PVTdataList.Pmark)./PVTeval3.PVTdataList.PErrTot)
PVTeval1 = fitErrModPVTeval(PVTeval1,[]);

opt = [];

% Initialize and perform fit
PVTeval1 = initPVTeval(name1,PVTdata,eosMgPvPrior1,opt);
PVTeval1 = fitSampModPVTeval(PVTeval1,[],fixFlag1);
% Review residuals to check if error model is needed
hist((PVTeval1.Psamp-PVTeval1.PVTdataList.Pmark)./PVTeval1.PVTdataList.PErrTot,15)
std((PVTeval1.Psamp-PVTeval1.PVTdataList.Pmark)./PVTeval1.PVTdataList.PErrTot)
PVTeval1 = fitErrModPVTeval(PVTeval1,[]);
hist((PVTeval1.Psamp-PVTeval1.PVTdataList.Pmark)./PVTeval1.PVTdataList.PErrTot,15)
std((PVTeval1.Psamp-PVTeval1.PVTdataList.Pmark)./PVTeval1.PVTdataList.PErrTot)
% Refit with updated error model
PVTeval1 = fitSampModPVTeval(PVTeval1,[],fixFlag1);
fileNm = '../data/Tange2012ErrTot.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval1.PVTdataList,fileNm);
fileNm = '../data/Tange2012Fit.md';
eosModList = [PVTeval1.sampEosFit,PVTeval1.sampEosPrior];
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'decimal');
%[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm,'paren');

% Initialize and perform fit
opt.robustFit = true;
PVTeval2 = initPVTeval(name2,PVTdata,eosMgPvPrior2,opt);
PVTeval2 = fitSampModPVTeval(PVTeval2,[],fixFlag2);
% Review residuals to check if error model is needed
hist((PVTeval2.Psamp-PVTeval2.PVTdataList.Pmark)./PVTeval1.PVTdataList.PErrTot,15)
std((PVTeval2.Psamp-PVTeval2.PVTdataList.Pmark)./PVTeval2.PVTdataList.PErrTot)
PVTeval2 = fitErrModPVTeval(PVTeval2,[]);
hist((PVTeval2.Psamp-PVTeval2.PVTdataList.Pmark)./PVTeval2.PVTdataList.PErrTot,15)
std((PVTeval2.Psamp-PVTeval2.PVTdataList.Pmark)./PVTeval2.PVTdataList.PErrTot)
% Refit with updated error model
PVTeval2 = fitSampModPVTeval(PVTeval2,[],fixFlag2);
PVTeval2.PVTdata
fileNm = '../data/Tange2012ErrTotR.md';
[tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTeval2.PVTdataList,fileNm);

fileNm = '../data/Tange2012Fit.md';
eosModList = [PVTeval2.sampEosFit,PVTeval2.sampEosPrior]
eosModLbl = {'fit','prior'};
[tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, eosModLbl,fileNm);


% Has std near 1

% fitErrMod and refit
PVTeval1 = fitErrModPVTeval(PVTeval1,errModInit);

% Test4: Restricted Thermal Fit No prior
name4 = 'Tange2012:Refine Test4 (restricted thermal fit, No prior)';
fixFlag4 = [1 1 1 1 0 0 1];
pEosPrior      = eosMgPvT12.pEos;
eosMgPvPrior4 = eosMgPvT12;
eosMgPvPrior4.pEos = pEosPrior;
eosMgPvPrior4.pEosCov = zeros(7);
eosMgPvPrior4.pEosCov(5,5) =   Inf;
eosMgPvPrior4.pEosCov(6,6) =   Inf;


PVTeval1 = initPVTeval(name1,PVTdata,eosMgPvPrior1,opt);
PVTeval2 = initPVTeval(name2,PVTdata,eosMgPvPrior2,opt);
PVTeval3 = initPVTeval(name3,PVTdata,eosMgPvPrior3,opt);
PVTeval4 = initPVTeval(name4,PVTdata,eosMgPvPrior4,opt);

%keyboard;
PVTeval1 = fitSampModPVTeval(PVTeval1,[],fixFlag1);
PVTeval2 = fitSampModPVTeval(PVTeval2,[],fixFlag2);
PVTeval3 = fitSampModPVTeval(PVTeval3,[],fixFlag3);
PVTeval4 = fitSampModPVTeval(PVTeval4,[],fixFlag4);


% Now test how reasonable the optimal estimation approximation is
%  Do this for Test2, since it fits all params except V0

pavg    = PVTeval2.sampEosFit.pEos;
pcov    = PVTeval2.sampEosFit.pEosCov;
fixFlag = PVTeval2.fixFlag;
costFun = PVTeval2.nLogPFun;

Ndraw = 3000;
[costDraw,costDrawApp,pdraw,costDrawMin] = ...
    testCostFun(Ndraw,pavg,pcov,fixFlag,costFun);
%logCostRatio = log((costDrawApp-costDrawMin)./(costDraw-costDrawMin));
costDiff = costDraw-costDrawApp;
wt = exp(-costDiff);

[costDrawS,indS] = sort(costDraw);
wtS = wt(indS);
cumWt = cumsum(wtS);
costDrawAppS = costDrawApp(indS);



plot(costDraw-costDrawMin,costDiff,'ko')
xlabel('CostFun - CostFunMin');
ylabel('CostFun Approx Error');

credLvl = [0.05:0.05:1];
inBoundFrac = zeros(size(credLvl));
outBoundCnt = zeros(size(credLvl));
for(i=1:length(credLvl))
    ilvl = credLvl(i);
    threshApp = interp1([1:Ndraw]/Ndraw,sort(costDrawApp-costDrawMin),ilvl);
    %thresh = interp1(cumWt/Ndraw,sort(costDrawS-costDrawMin),ilvl)
    indCredLvl = find(cumWt/cumWt(end)<= ilvl);
    inBoundFrac(i) = mean(costDrawAppS(indCredLvl)-costDrawMin<=threshApp);
    outBoundCnt(i) = sum(costDrawAppS(indCredLvl)-costDrawMin<=threshApp)/Ndraw;
end
plot(100*credLvl,-100*(outBoundCnt-credLvl),'ro-');
xlabel('Confidence Level [%]','FontSize',16);
ylabel('Error in Confidence Level [%]','FontSize',16)


plot(costDrawAppS-costDrawMin, [1:Ndraw]/Ndraw,'k-',...
    costDrawS-costDrawMin,cumWt/Ndraw,'r-')
xlabel('CostFun - CostFunMin');
ylabel('CDF');
legend('Norm-Approx','Reweighted','Location','SouthEast');

plot(costDrawS-costDrawMin,cumWt/Ndraw,'r-',...
    sort(costDrawApp)-costDrawMin, [1:Ndraw]/Ndraw,'k-')




[pFree,pFree,pcovFree] = getFreeParams(p,p,pcov,fixFlag);
[pAll] = getAllParams(pFree,p,fixFlag);

sqrt(diag(PVTeval2.sampEosFit.pEosCov))'

PVTeval1.nLogPFun(PVTeval1.sampEosFit.pEos)

% costFun compare

