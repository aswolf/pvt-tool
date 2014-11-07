% setPVTdataStd - load PVTdata struc using (P,V,T) with their errors
%
%  'std' error mode assumes that errors on P,V,T are independent
%       * This is not actually true, due to correlations between T and P
%       * errors are added in quadrature estimated from current samp Eos model
function PVTdata = setPVTdataStd(PVTdata,P,V,T,PErr,VErr,TErr,measGrpID)

    errMode = 'std';

    if(isempty(measGrpID))
        measGrpID = {'1'};
    end
    Ndat = length(V);
    assert(size(P,1) == Ndat,'P must have one row for each datum')
    assert(size(T,1) == Ndat,'T must have one entry for each datum')
    assert(size(PErr,1) == Ndat,'PErr must have one row for each datum')
    assert(size(VErr,1) == Ndat,'VErr must have one entry for each datum')
    assert(size(TErr,1) == Ndat,'TErr must have one entry for each datum')
    assert(iscellstr(measGrpID),'measGrpID must be a cell array of str')

    PVTdata.errMode  = errMode;

    PVTdata.Pmark = P;
    PVTdata.V     = V;
    PVTdata.T     = T;
    PVTdata.PErr  = PErr;
    PVTdata.VErr  = VErr;
    PVTdata.TErr  = TErr;

    PVTdata.measGrpID= measGrpID;   
end
