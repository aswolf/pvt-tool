% setPVTdataMark - load PVTdata struc using (Vmark,V,T) with their errors
%
%  'mark' error mode treats errors on Vmark,V,T as independent
%       * This is realistic!
%       * errors are added in quadrature estimated from mark and samp Eos
%       models
function PVTdata = setPVTdataMark(PVTdata,markLbl,markEos,...
        Vmark,V,T,VmarkErr,VErr,TErr,measGrpID)

    errMode = 'mark';

    if(isempty(measGrpID))
        measGrpID = {'1'};
    end
    Ndat = length(V);
    assert(size(Vmark,1) == Ndat,'Vmark must have one row for each datum')
    assert(size(T,1) == Ndat,'T must have one entry for each datum')
    assert(size(VmarkErr,1) == Ndat,'VmarkErr must have one row for each datum')
    assert(size(VErr,1) == Ndat,'VErr must have one entry for each datum')
    assert(size(TErr,1) == Ndat,'TErr must have one entry for each datum')
    assert(iscellstr(measGrpID),'measGrpID must be a cell array of str')


    PVTdata.errMode  = errMode;
    PVTdata.markLbl = markLbl;
    PVTdata.markEos = markEos;

    [Pmark,PmarkDerivs] = evalPressEos([],markEos,Vmark,T);

    PVTdata.Pmark = Pmark;
    PVTdata.Vmark = Vmark;
    PVTdata.V     = V;
    PVTdata.T     = T;
    PVTdata.VmarkErr = VmarkErr;
    PVTdata.VErr     = VErr;
    PVTdata.TErr     = TErr;

    PVTdata.measGrpID= measGrpID;   

    PVTdata.PmarkDerivs = PmarkDerivs;
end
