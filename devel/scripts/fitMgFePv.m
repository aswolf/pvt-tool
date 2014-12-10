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



