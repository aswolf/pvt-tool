function PVTeval = fitSampModPVTeval(PVTeval)
    % fit sample eos given total errors defined in PVTdata
    % update sampEos in PVTeval
    % updateModPVTeval (with new presid, sampPderivs, Perr
    %    - should this be a separate function?
    %

    eosMod = PVTeval.sampEos;

    T0                 = eosMod.T0                 ;
    NpCold             = eosMod.NpCold             ;
    coldEosFun         = eosMod.coldEosFun         ;
    hotEosFun          = eosMod.hotEosFun          ;
    hotExtraInputs     = eosMod.hotExtraInputs     ;
    addedThermPressFun = eosMod.addedThermPressFun ;

    priorEos = eosMod.pEos;
    pinitEos = priorEos;


    [pfit pfitcov nLogPFun] = fitHotCompressData(pinitEos,fixFlag,...
        T0,NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun,P,V,T,PerrTot,opt);
    
    pressFun = @(V,T,pEos)(calcPressThermAddEos(V,T,T0,...
        pEos(1:NpCold),pEos(NpCold+1:end),coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun));




    pressFun = @(V,T,pEos)(calcPressThermAddEos(V,T,T0,...
        pEos(1:NpCold),pEos(NpCold+1:end),coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun));
    [P,KT,Cv,gam,thmExp] = pressFun(V,T,pEos);
end
