%
%
%    PVTdata.PResid       = [];
%    PVTdata.PSampDerivs  = [];
%    PVTdata.PMarkDerivs  = [];
%
    % Determine initial errMod and errModCov
    %[measGrpID,errMod,errModCov] = initErrMod(measGrpID,Ndat);
    %if(length(measGrpID)==1)
    %    val = measGrpID;
    %    measGrpID = val*ones(size(P));
    %end
function PVTdata = initPVTdata(name,Pmark,V,T,PmarkErr,VErr,TErr,measGrpID,...
        markLbl,markEos,Vmark,VmarkErr,opt)

    optDefault = getPVTdataDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    assert(length(markEos) <= 1,'Multiple Press markers NOT YET implimented');
    assert(strcmp(class(markLbl),'char') | ...
        (strcmp(class({markLbl}),'cell')& length(markLbl) <= 1),...
        'Multiple Press markers NOT YET implimented');
    assert(isstr(name),'name must be a string name for this PVT dataset.')

    if(isempty(measGrpID))
        measGrpID = 1;
    end

    calcPmark = true;
    invVmark  = false;
    pressFixed = false;
    if(isempty(Pmark))
        assert(~isempty(markEos),...
            'If Pmark is empty, must provide markEos to calculate press.');
    end
    if(isempty(Vmark))
        assert(~isempty(markEos),...
            'If Vmark is empty, must provide markEos to invert for Vmark.');
        invVmark = true;
    end
    if(isempty(markEos))
        calcPmark = false;
        pressFixed = true;
    end
    assert(invVmark==false,['Invert Vmark option is NOT yet implimented. '...
        'Give Vmark explicitly or impliment this option in initPVTdata.']);

    Ndat = length(V);

    PVTdata.name     = name;

    PVTdata.markLbl  = markLbl;
    PVTdata.markEos  = markEos;
    PVTdata.pressFixed= pressFixed;

    PVTdata.measGrpID= measGrpID;   

    PVTdata.Pmark    = Pmark;   
    PVTdata.V        = V;  
    PVTdata.T        = T;  
    PVTdata.PmarkErr = PmarkErr;
    PVTdata.VErr     = VErr;
    PVTdata.TErr     = TErr;

    PVTdata.Vmark    = Vmark; 
    PVTdata.VmarkErr = VmarkErr;    
    PVTdata.Pderivs = [];    

    if(calcPmark)
        PTOLABS = opt.PTOLABS;
        [Pcalc,Pderivs] = evalPressEos([],markEos,Vmark,T);
        assert(all(abs(Pcalc-Pmark)<PTOLABS), ...
            'Pressure from markEos and input must agree (within TOL in opt)');
        PVTdata.Pmark = Pcalc;
        PVTdata.Pderivs = Pderivs;
    end
end
function [measGrpID,errMod,errModCov] = initErrMod(measGrpID,Ndat,opt)
    uniqMeasGrpID = unique(measGrpID);
    if(length(measGrpID)==1)
        val = measGrpID;
        measGrpID = val*ones(Ndat,1);
    end
    errMod = ones(length(uniqMeasGrpID),3);
    errModPrior = [opt.errModVmarkPrior opt.errModVmarkPrior opt.errModVmarkPrior];
    errModCov = diag(Inf*ones(1,numel(errMod)));
end
