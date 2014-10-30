function PVTdata = updatePVTdata(PVTdata)
    markEos = PVTdata.markEos;
    Vmark   = PVTdata.Vmark;
    T       = PVTdata.T;

    assert(~isempty(markEos),'Eos Press marker must be defined.');
    assert(length(markEos) == 1,'Multiple Press markers NOT YET implimented');
    [Pcalc,Pderivs] = evalPressEos([],markEos,Vmark,T);

    PVTdata.Pmark = Pcalc;
    PVTdata.Pderivs = Pderivs;
end
