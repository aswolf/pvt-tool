% readPVTdataTable - standard method for reading PVT data files
function PVTdata = readPVTdataTable(fileName,name,measGrpIDLbl,errMode,...
        markLbl,markEos,opt)

    optDefault = getPVTdataDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    errModeList ={'mark','std','tot'}; 
    assert( any(strcmp(errModeList,errMode)), ...
        ['errMode must be from recognized options: ' strjoin(errModeList)])
    nheadlines = 3;
    dataRead = importdata(fileName,'|',nheadlines);
    dataTbl = dataRead.data;

    % Column definition depends on errMode
    switch errMode
        case 'mark'
            % measGrpInd T TErr Vmark VmarkErr V VErr
            measGrpInd = dataTbl(:,1);
            T          = dataTbl(:,2);
            TErr       = dataTbl(:,3);
            Vmark      = dataTbl(:,4);
            VmarkErr   = dataTbl(:,5);
            V          = dataTbl(:,6);
            VErr       = dataTbl(:,7);
        case 'std'
            % measGrpInd T TErr P PErr V VErr
            measGrpInd = dataTbl(:,1);
            T          = dataTbl(:,2);
            TErr       = dataTbl(:,3);
            P          = dataTbl(:,4);
            PErr       = dataTbl(:,5);
            V          = dataTbl(:,6);
            VErr       = dataTbl(:,7);
        case 'tot'
            % measGrpInd T P V PErrTot
            measGrpInd = dataTbl(:,1);
            T          = dataTbl(:,2);
            P          = dataTbl(:,3);
            PErrTot    = dataTbl(:,4);
            V          = dataTbl(:,5);
    end


    % Scale volumes if option selected
    if(~isnan(opt.VScl))
        V    = opt.VScl*V;
        VErr = opt.VScl*VErr;
    end
    if(~isnan(opt.VmarkScl))
        Vmark    = opt.VmarkScl*Vmark;
        VmarkErr = opt.VmarkScl*VmarkErr;
    end

    uniqMeasGrpInd = unique(measGrpInd);
    assert(all(uniqMeasGrpInd >= 0 & mod(uniqMeasGrpInd,1)==0),...
        'measGrpInd must all be non-negative integers');
    assert(max(uniqMeasGrpInd) <= length(measGrpIDLbl),...
        'measGrpIDLbl must provide an ID label for every measGrpInd')

    % mask out zero-values for meas grp Ind 
    mask = measGrpInd==0;

    measGrpID = cell([length(V),1]);
    for(i=1:length(V))
        if(mask(i))
            measGrpID(i) = {NaN};
        else
            measGrpID(i) = measGrpIDLbl(measGrpInd(i));
        end
    end

    PVTdata = initPVTdata(name,opt);

    switch errMode
        case 'mark'
            Vmark = Vmark(~mask);
            V = V(~mask);
            T = T(~mask);
            VmarkErr = VmarkErr(~mask);
            VErr = VErr(~mask);
            TErr = TErr(~mask);
            measGrpID = measGrpID(~mask);
            PVTdata = setPVTdataMark(PVTdata,markLbl,markEos,...
                Vmark,V,T,VmarkErr,VErr,TErr,measGrpID);
        case 'std'
            P = P(~mask);
            V = V(~mask);
            T = T(~mask);
            PErr = PErr(~mask);
            VErr = VErr(~mask);
            TErr = TErr(~mask);
            measGrpID = measGrpID(~mask);
            PVTdata = setPVTdataStd(PVTdata,P,V,T,PErr,VErr,TErr,measGrpID);
        case 'tot'
            P = P(~mask);
            V = V(~mask);
            T = T(~mask);
            PErrTot = PErrTot(~mask);
            measGrpID = measGrpID(~mask);
            %disp('NOTE: this option is just for proof of concept! Do not use')
            PVTdata = setPVTdataTot(PVTdata,P,V,T,PErrTot,measGrpID);
    end

end

