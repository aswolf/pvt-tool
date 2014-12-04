% loadPVTdata - function to load PVTdata from loadScript
function PVTdata = loadPVTdata(loadScriptNm,markEos)
    [dataDir,fileName,name,measGrpIDLbl,errMode,markLbl,optDefault] = ...
        initLoadVars(loadScriptNm);
    opt = optDefault;

    [path,file] = fileparts(loadScriptNm);

    eval(file);
    opt = setDefaultOpt(opt,optDefault);
    fileName = fullfile(dataDir,fileName);

    PVTdata = readPVTdataTable(fileName,name,measGrpIDLbl,errMode,...
        markLbl,markEos,opt);
end
function [dataDir,fileName,name,measGrpIDLbl,errMode,markLbl,optDefault] = ...
        initLoadVars(loadScriptNm)
    dataDir = fileparts(loadScriptNm);
    fileName = [];
    name = [];
    measGrpIDLbl = {};
    errMode = [];
    markLbl = '';
    optDefault = getPVTdataDefaultOpt();
end
