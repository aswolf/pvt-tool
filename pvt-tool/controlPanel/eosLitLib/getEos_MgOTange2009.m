function eosMod = getEos_MgOTange2009()
    eosMod = getBestModel();
end
function eosMod = getBestModel()
    material = 'MgO';
    name = 'MgO-Tange2009';
    T0 = 300;
    V0 = 74.698;
    K0 = 160.63;
    KP0= 4.367;

    Natom = 4*2;
    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 a b 1.0];

    pEos = [pColdEos pHotEos];
    NpCold = length(pColdEos);

    % DOES NOT INCLUDE CORRELATION
    pEosCov = diag([0 0.18 0.013 13 0.015 0.019 1.1 0]);

    coldEosFun = @VinetEos;
    hotEosFun = @MieGrunDebyeHotEos;
    debyeDerivsFun = @debyeTange;
    
    hotExtraInputs = {Natom, debyeDerivsFun};
    addedThermPressFun = [];

    eosMod = initEos(name,material,pEos,pEosCov,T0,NpCold,coldEosFun,hotEosFun,hotExtraInputs,addedThermPressFun);
end
