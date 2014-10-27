function eosMod = initEos(pEos,pEosCov,T0,NpCold,coldEosFun,hotEosFun,...
        hotExtraInputs,addedThermPressFun)

    eosMod.pEos               = pEos;
    eosMod.pEosCov            = pEosCov;
    eosMod.T0                 = T0;
    eosMod.NpCold             = NpCold;
    eosMod.coldEosFun         = coldEosFun;
    eosMod.hotEosFun          = hotEosFun;
    eosMod.hotExtraInputs     = hotExtraInputs;
    eosMod.addedThermPressFun = addedThermPressFun;
end
