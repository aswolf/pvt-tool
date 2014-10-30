function [P,Pderivs,KT,Cv,gam,thmExp] = evalPressEos(pEos,eosMod,V,T)
    if(isempty(pEos))
        pEos = eosMod.pEos;
    end
    assert(size(V,2)==1,'V must be a col vector');
    assert(size(T,2)==1,'T must be a col vector');

    T0                 = eosMod.T0                 ;
    NpCold             = eosMod.NpCold             ;
    coldEosFun         = eosMod.coldEosFun         ;
    hotEosFun          = eosMod.hotEosFun          ;
    hotExtraInputs     = eosMod.hotExtraInputs     ;
    addedThermPressFun = eosMod.addedThermPressFun ;

    pressFun = @(V,T,pEos)(calcPressThermAddEos(V,T,T0,...
        pEos(1:NpCold),pEos(NpCold+1:end),coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun));
    [P,KT,Cv,gam,thmExp] = pressFun(V,T,pEos);
    dPdV = -KT./V;
    dPdT = KT.*thmExp;
    Pderivs = [dPdV,dPdT];
end
