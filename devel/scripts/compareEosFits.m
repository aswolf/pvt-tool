% Perform series of fits under different conditions to compare with MINUTI
% results

% Use Tange 2012 pv dataset
%PVTdata = getPVTdata_MgPvTange2012();
PVTdata = getPVTdata_MgPvTange2012_compare();
% load Tange Eos for MgPv and MgO
eosMgPvT12 = getEos_MgPvTange2012();
eosMgOMod = getEos_MgOTange2009();

% Perform multiple tests with different priors

% Test1: Prior all Tange values
name1 = 'Tange2012:Refine Test1 (Tange prior)';
fixFlag1 = [1 0 0 0 0 0 1];
pEosPrior      = eosMgPvT12.pEos;
eosMgPvPrior1 = eosMgPvT12;
eosMgPvPrior1.pEos = pEosPrior;
eosMgPvPrior1.pEosCov = zeros(7);
eosMgPvPrior1.pEosCov(2,2) = Inf;
eosMgPvPrior1.pEosCov(3,3) = Inf;
eosMgPvPrior1.pEosCov(4,4) = 200^2;
eosMgPvPrior1.pEosCov(5,5) =   1^2;
eosMgPvPrior1.pEosCov(6,6) =   1^2;

% Test2: Weakly Informed Thermal Prior
name2 = 'Tange2012:Refine Test2 (Weakly informed thermal prior)';
fixFlag2 = [1 0 0 0 0 0 1];
pEosPrior      = eosMgPvT12.pEos;
pEosPrior(4) = 1000;
pEosPrior(5:6) = 1;
eosMgPvPrior2 = eosMgPvT12;
eosMgPvPrior2.pEos = pEosPrior;
eosMgPvPrior2.pEosCov = zeros(7);
eosMgPvPrior2.pEosCov(2,2) = Inf;
eosMgPvPrior2.pEosCov(3,3) = Inf;
eosMgPvPrior2.pEosCov(4,4) = 200^2;
eosMgPvPrior2.pEosCov(5,5) =   1^2;
eosMgPvPrior2.pEosCov(6,6) =   1^2;

% Test3: Restricted Thermal Fit
name3 = 'Tange2012:Refine Test3 (restricted thermal fit)';
fixFlag3 = [1 1 1 1 0 0 1];
pEosPrior      = eosMgPvT12.pEos;
pEosPrior(5:6) = 1;
eosMgPvPrior3 = eosMgPvT12;
eosMgPvPrior3.pEos = pEosPrior;
eosMgPvPrior3.pEosCov = zeros(7);
eosMgPvPrior3.pEosCov(5,5) =   1^2;
eosMgPvPrior3.pEosCov(6,6) =   1^2;

% Test4: Restricted Thermal Fit No prior
name4 = 'Tange2012:Refine Test3 (restricted thermal fit, No prior)';
fixFlag4 = [1 1 1 1 0 0 1];
pEosPrior      = eosMgPvT12.pEos;
eosMgPvPrior4 = eosMgPvT12;
eosMgPvPrior4.pEos = pEosPrior;
eosMgPvPrior4.pEosCov = zeros(7);
eosMgPvPrior4.pEosCov(5,5) =   Inf;
eosMgPvPrior4.pEosCov(6,6) =   Inf;

opt = [];

PVTeval1 = initPVTeval(name1,PVTdata,eosMgPvPrior1,opt);
PVTeval2 = initPVTeval(name2,PVTdata,eosMgPvPrior2,opt);
PVTeval3 = initPVTeval(name3,PVTdata,eosMgPvPrior3,opt);
PVTeval4 = initPVTeval(name4,PVTdata,eosMgPvPrior4,opt);

%keyboard;
PVTeval1 = fitSampModPVTeval(PVTeval1,[],fixFlag1);
PVTeval2 = fitSampModPVTeval(PVTeval2,[],fixFlag2);
PVTeval3 = fitSampModPVTeval(PVTeval3,[],fixFlag3);
PVTeval4 = fitSampModPVTeval(PVTeval4,[],fixFlag4);

PVTeval1.

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

