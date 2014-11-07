% updatePVTdata - update PErrTot appropriate to selected errorMode
function PVTdata = updatePVTdata(PVTdata,PsampDerivs)
    errMode  = PVTdata.errMode;

    switch(errMode)
        case 'std'
            assert(~isempty(PsampDerivs),'PsampDerivs must be given for std err mode.');
            P  = PVTdata.Pmark;
            V  = PVTdata.V;
            T  = PVTdata.T;

            PErr  = PVTdata.PErr;
            VErr  = PVTdata.VErr;
            TErr  = PVTdata.TErr;

            PErrTot = sqrt(PErr.^2 + (VErr.*PsampDerivs(:,1)).^2 + ...
                (TErr.*PsampDerivs(:,2)).^2);

            PVTdata.PErrTot = PErrTot;
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


            PErrTot = sqrt((VmarkErr.*PmarkDerivs(:,1)).^2 + ...
                (VErr.*PsampDerivs(:,1)).^2 + ...
                (TErr.*(PmarkDerivs(:,2)-PsampDerivs(:,2))).^2);

            PVTdata.PErrTot = PErrTot;
        case 'tot'
            % No need to do anything since total error is fixed
        otherwise
            assert(false,[errMode 'is not a recognized error Mode option']);
    end
end
