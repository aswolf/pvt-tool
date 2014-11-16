function PVTdata = getPVTdata_MgFePvWolf2014()

    errMode = 'mark';
    neline = '111';

    %errModeList ={'mark','std','tot'}; 
    %if(~exist('errMode'))
    %    errMode = 'mark';
    %end
    %assert( any(strcmp(errModeList,errMode)), ...
    %    ['errMode must be from recognized options: ' strjoin(errModeList)])
    %   1   2,     3      4, 5           6, 7      8,     9 
    %   grp P-VIN, err    V, err         T, err    V-MgO, err  
    
    VVTdat = importdata('Wolf2014MgFePvrawdat.txt');
    VVT = VVTdat.data;
    %   1   2,     3      4, 5           6, 7      8,     9 
    %   grp P-VIN, err    V, err         T, err    V-MgO, err  
    V     = VVT(:,1);
    VErr  = VVT(:,2);
    T     = VVT(:,3);
    TErr  = VVT(:,4);
    Vne111= VVT(:,5);
    Vne200= VVT(:,6);

    if(strcmp(neline,'111'))
        Vmark = Vne111;
    elseif(strcmp(neline,'200'))
        Vmark = Vne200;
    end

    % Assume that typical fractional error on volume
    % NOTE: Major assumption.
    % Initial guess for marker error volume is that 
    VErrFac = exp(median(log(VErr./V)));
    VmarkErr = VErrFac*Vmark;

    measGrpInd = ones(size(V));
    measGrpInd(T>310) = 2;

    measGrpID = cellstr(num2str(measGrpInd));

    % Tange 2012 uses MgO press scl of Tange 2009
    name = 'MgFePvData-Wolf2014';

    eosMod_Ne = getEos_NeDewaele2008();

    markLbl = eosMod_Ne.name;
    opt = [];

    PVTdata = initPVTdata(name,opt);

    switch errMode
        case 'mark'
            PVTdata = setPVTdataMark(PVTdata,markLbl,eosMod_Ne,...
                Vmark,V,T,VmarkErr,VErr,TErr,measGrpID);
        case 'std'
            PVTdata = setPVTdataStd(PVTdata,P,V,T,PErr,VErr,TErr,measGrpID);
        case 'tot'
            disp('NOTE: this option is just for proof of concept! Do not use')
            PerrTot = sqrt(PErr.^2 + 1);
            PVTdata = setPVTdataTot(PVTdata,P,V,T,PErrTot,measGrpID);
    end
end

