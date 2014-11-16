function eosMod = getEos_NeDewaele2008()
    eosMod = getBestModel();
end
function eosMod = getBestModel()
    material = 'Ne';
    name = 'Ne-Dewaele2008';

    Natom = 4;

    % NOTE: Dewaele's EOS has an absolute zero reference temp
    T0 = 0;
    V0 = Natom*22.234;
    K0 = 1.070;
    KP0= 8.40;

    Tdeb0 = 75.1;
    gam1  = 2.442;
    q = 0.97;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam1 q 1.0];
    
    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    % DOES NOT INCLUDE CORRELATION
    pEosCov = diag([0 0.016 0.03 0 0 0 0]);

    coldEosFun = @VinetEos;
    hotEosFun = @MieGrunDebyeHotEos;
    debyeDerivsFun = @debyeShiftPowerLaw;

    hotExtraInputs = {Natom, debyeDerivsFun};
    addedThermPressFun = [];

    eosMod = initEos(name,material,pEos,pEosCov,T0,...
        NpCold,coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);
end
