% setPVTdataTot - load PVTdata struc using (P,V,T) with total P error estimate
%
%  'tot' error mode enables the user to set and fix PErrTot manually
function PVTdata = setPVTdataTot(PVTdata,P,V,T,PErrTot,measGrpID)

    errMode = 'tot';

    if(isempty(measGrpID))
        measGrpID = {'1'};
    end
    Ndat = length(V);
    assert(size(P,1) == Ndat,'P must have one row for each datum')
    assert(size(T,1) == Ndat,'T must have one entry for each datum')
    assert(size(PErrTot,1) == Ndat,'PErrTot must have one row for each datum')
    assert(iscellstr(measGrpID),'measGrpID must be a cell array of str')

    PVTdata.errMode  = errMode;

    PVTdata.Pmark = P;
    PVTdata.V     = V;
    PVTdata.T     = T;
    PVTdata.PErrTot = PErrTot;

    PVTdata.measGrpID= measGrpID;   

    PVTdata.PErrTot  = PErrTot;
end
