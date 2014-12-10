dataDir = '/Users/aswolf/Documents/code/MATLAB/pvt-tool/devel/data';

% Fit MgPv (Tange2012) and MgFePv (Wolf2014) datasets
lsqMgPv  =runPVTtool(fullfile(dataDir,'fit_Tange2012_lsq_fixdeb.in'));
lsqMgFePv=runPVTtool(fullfile(dataDir,'fit_Wolf2014_lsq_fixdeb.in'));

viewPVTFit(lsqMgFePv,'normal')
viewPVTFit(lsqMgFePv,'reduced')
viewPVTFit(lsqMgFePv,'hist')

viewPVTFit(lsqMgPv,'normal')
viewPVTFit(lsqMgPv,'reduced')
viewPVTFit(lsqMgPv,'hist')


Pbnds = [30, 130];
Tbnds = [300, 3000];
eosModEndmem = [lsqMgPv.sampEosFit, lsqMgFePv.sampEosFit];
Xendmem = [0, 0.13];
Xideal = [0:.02:0.2];
eosModIdeal = calcIdealLatticeMix(Xideal, Xendmem, eosModEndmem,Pbnds,Tbnds);

pEosIdealFit = reshape([eosModIdeal.pEos],[],length(eosModIdeal))';
plot(Xideal,pEosIdealFit(:,1),'k-')
plot(Xideal,pEosIdealFit(:,2),'k-')
plot(Xideal,pEosIdealFit(:,3),'k-')
plot(Xideal,pEosIdealFit(:,5),'k-')
plot(Xideal,pEosIdealFit(:,6),'k-')
