function eosMod = initEos(name,material,pEos,pEosCov,T0,NpCold,coldEosFun,hotEosFun,...
        hotExtraInputs,addedThermPressFun)
    assert(isstr(name),'name must be a string');
    assert(isstr(material),'material must be a string');

    eosMod.name               = name;
    eosMod.material           = material;
    eosMod.pEos               = pEos;
    eosMod.pEosCov            = pEosCov;
    eosMod.T0                 = T0;
    eosMod.NpCold             = NpCold;
    eosMod.coldEosFun         = coldEosFun;
    eosMod.hotEosFun          = hotEosFun;
    eosMod.hotExtraInputs     = hotExtraInputs;
    eosMod.addedThermPressFun = addedThermPressFun;
end
