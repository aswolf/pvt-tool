function [jointDatStruc] = joinEOSDataSets(datStruc)
    Ndat = length(datStruc);
    Nmeas = length(cat(1,datStruc.P));

    P        = NaN*ones(Nmeas,1);
    V        = NaN*ones(Nmeas,1);
    T        = NaN*ones(Nmeas,1);
    PerrTot  = NaN*ones(Nmeas,1);
    datsetID = NaN*ones(Nmeas,1);
    runID    = NaN*ones(Nmeas,1);
    residP   = NaN*ones(Nmeas,1);

    datsetName = {datStruc.datsetName};

    ind=1;
    for(i = 1:Ndat)
        idatStruc = datStruc(i);
        iNmeas = length(idatStruc.P);

        datsetID(ind:ind+iNmeas-1) = i;
        runID(ind:ind+iNmeas-1) = idatStruc.runID;
        P(ind:ind+iNmeas-1) = idatStruc.P;
        V(ind:ind+iNmeas-1) = idatStruc.V;
        T(ind:ind+iNmeas-1) = idatStruc.T;
        PerrTot(ind:ind+iNmeas-1) = idatStruc.PerrTot;
        ind = ind + iNmeas;
    end

    jointDatStruc.datsetName = datsetName;
    jointDatStruc.datsetID = datsetID;
    jointDatStruc.runID = runID;
    jointDatStruc.P = P;
    jointDatStruc.V = V;
    jointDatStruc.T = T;
    jointDatStruc.PerrTot = PerrTot;
end
