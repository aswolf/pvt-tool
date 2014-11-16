function PVTeval = fitSampModPVTeval(PVTeval,eosInitMod,fixFlag)
    % fit sample eos given total errors defined in PVTdata
    % update sampEos in PVTeval
    %
    eosPriorMod = PVTeval.sampEosPrior;
    if(isempty(eosInitMod))
        eosInitMod = eosPriorMod;
    end

    PVTdataList = PVTeval.PVTdataList;
    assert(numel(PVTdataList)==1,...
        'Multiple datasets in PVTdataList not yet implemented');
    PVTdata = PVTdataList(1);

    P = PVTdata.Pmark;
    V = PVTdata.V;
    T = PVTdata.T;

    % Propagate errors to total errors in pressure
    errModList = PVTeval.errModList;
    [Psamp,PsampDerivs] = evalPressEos([],eosInitMod,V,T);
    PVTdata = updatePVTdata(PVTdata,PsampDerivs,errModList);
    PerrTot = PVTdata.PErrTot;

    T0                 = eosPriorMod.T0                 ;
    NpCold             = eosPriorMod.NpCold             ;
    coldEosFun         = eosPriorMod.coldEosFun         ;
    hotEosFun          = eosPriorMod.hotEosFun          ;
    hotExtraInputs     = eosPriorMod.hotExtraInputs     ;
    addedThermPressFun = eosPriorMod.addedThermPressFun ;

    priorEos = eosPriorMod.pEos;
    pinitEos = eosInitMod.pEos;
    priorcovEos = eosPriorMod.pEosCov;
    pinitcovEos = eosPriorMod.pEosCov;

    opt = PVTeval.opt;

    [pfitEos pfitcovEos nLogPFun PressTotFun opt] = fitHotCompressData(pinitEos,fixFlag,...
        T0,NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun,P,V,T,PerrTot,opt);

    eosFitMod = eosInitMod;
    eosFitMod.pEos = pfitEos;
    eosFitMod.pEosCov = pfitcovEos;

    %pressFun = @(V,T,pEos)(calcPressThermAddEos(V,T,T0,...
    %    pEos(1:NpCold),pEos(NpCold+1:end),coldEosFun,hotEosFun,hotExtraInputs,...
    %    addedThermPressFun));
    [Psamp,KTsamp,Cvsamp,gamsamp,thmExpsamp] = PressTotFun(V,T,pfitEos);

    dPsampdV    = -KTsamp./V;
    dPsampdT    = KTsamp.*thmExpsamp;
    PsampDerivs = [dPsampdV,dPsampdT];

    % Propagate errors to total errors in pressure
    PVTdata = updatePVTdata(PVTdata,PsampDerivs,errModList);

    %Update list element
    PVTdataList(1) = PVTdata;

    PVTeval.PVTdataList  = PVTdataList;
    PVTeval.fixFlag      = fixFlag;
    PVTeval.sampEosInit = eosInitMod;
    PVTeval.sampEosFit   = eosFitMod;
    PVTeval.nLogPFun  = nLogPFun;
    PVTeval.opt = opt;
    PVTeval.Psamp = Psamp;
    PVTeval.PsampDerivs = PsampDerivs;
end
