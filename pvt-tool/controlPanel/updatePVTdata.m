% updatePVTdata - update PErrTot appropriate to selected errorMode
function PVTdata = updatePVTdata(PVTdata,PsampDerivs,errModList)
    errMode  = PVTdata.errMode;
    measGrpID = PVTdata.measGrpID;
    uniqID = unique(measGrpID);


    switch(errMode)
        case 'std'
            assert(~isempty(PsampDerivs),'PsampDerivs must be given for std err mode.');
            P  = PVTdata.Pmark;
            V  = PVTdata.V;
            T  = PVTdata.T;

            PErr  = PVTdata.PErr;
            VErr  = PVTdata.VErr;
            TErr  = PVTdata.TErr;

            VsampErrFac = zeros(size(VErr));
            TErrFac = zeros(size(TErr));
            for(i=1:length(uniqID))
                ierrMod = errModList(i,:);
                imeasGrpInd = find(strcmp(measGrpID,uniqID{i}));
                VsampErrFac(imeasGrpInd) = exp(ierrMod(1));
                TErrFac(imeasGrpInd)     = exp(ierrMod(2));
            end
            PErrTerms = zeros(length(V),3);
            PErrTerms(:,1) = PErr;
            PErrTerms(:,2) = abs(VsampErrFac.*VErr.*PsampDerivs(:,1));
            PErrTerms(:,3) = abs(TErrFac.*TErr.*PsampDerivs(:,2));
            PErrTot = sqrt(sum(PErrTerms.^2,2));

            PVTdata.PErrTot = PErrTot;
            PVTdata.PErrTerms = PErrTerms;
        case 'mark'
            assert(~isempty(PsampDerivs),'PsampDerivs must be given for std err mode.');

            markEos = PVTdata.markEos;
            Vmark = PVTdata.Vmark;
            V     = PVTdata.V;
            T     = PVTdata.T;

            VmarkErr = PVTdata.VmarkErr;
            VErr     = PVTdata.VErr;
            TErr     = PVTdata.TErr;

            PmarkDerivs = PVTdata.PmarkDerivs;



            VmarkErrFac = zeros(size(VmarkErr));
            VsampErrFac = zeros(size(VErr));
            TErrFac = zeros(size(TErr));
            for(i=1:length(uniqID))
                ierrMod = errModList(i,:);
                imeasGrpInd = find(strcmp(measGrpID,uniqID{i}));
                VmarkErrFac(imeasGrpInd) = exp(ierrMod(1));
                VsampErrFac(imeasGrpInd) = exp(ierrMod(1));
                TErrFac(imeasGrpInd)     = exp(ierrMod(2));
            end

            PErrTerms = zeros(length(V),3);
            PErrTerms(:,1) = abs(VmarkErrFac.*VmarkErr.*PmarkDerivs(:,1));
            PErrTerms(:,2) = abs(VsampErrFac.*VErr.*PsampDerivs(:,1));
            PErrTerms(:,3) = abs(TErrFac.*TErr.*(PmarkDerivs(:,2)-PsampDerivs(:,2)));
            PErrTot = sqrt(sum(PErrTerms.^2,2));


            PVTdata.PErrTot = PErrTot;
            PVTdata.PErrTerms = PErrTerms;
        case 'tot'
            % No need to do anything since total error is fixed
        otherwise
            assert(false,[errMode 'is not a recognized error Mode option']);
    end
end
