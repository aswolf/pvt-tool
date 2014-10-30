function eosMod = getEos_MgPvTange2012()
    eosMod = getBestModel();
end
function eosMod = getBestModel()
    material = 'MgPv';
    name = 'MgPv-Tange2012';
    T0 = 300;

    V0 = 162.373;
    K0 = 258.4;
    KP0= 4.10;

    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyePowerLaw;
    Natom = 4*5;
    Tdeb0 = 940;
    gam0  = 1.55;
    q = 1.1;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 q 1.0];
    
    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    % DOES NOT INCLUDE CORRELATION
    pEosCov = diag([0 1.7 0.07 140 0.09 0.3 0]);

    coldEosFun = @VinetEos;
    hotEosFun = @MieGrunDebyeHotEos;
    debyeDerivsFun = @debyePowerLaw;

    hotExtraInputs = {Natom, debyeDerivsFun};
    addedThermPressFun = [];

    eosMod = initEos(name,material,pEos,pEosCov,T0,NpCold,coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);
end
