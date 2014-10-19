function [thermExpTrend,VTrend,TTrend] = calcThermExpTrend(Ptarget,...
        Tmin,Tmax,Nstep,T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,hotExtraInputs,elecThermPressFun)

    assert(Tmin == T0,'Code not generalized to non-ambient starting Temp.');

    V0 = pColdEos(1);
    Vnext = V0;
    TTrend = linspace(Tmin,Tmax,Nstep);
    dT = diff(TTrend(1:2));
    PTrend = zeros(1,Nstep);
    for(i=1:Nstep)
        Tcurr = TTrend(i);
        Vcurr = Vnext;
        VTrend(i) = Vcurr;
        [iP,iKT,iCv,igam] = calcPressThermAddEos(Vcurr,Tcurr,T0,...
            pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
            elecThermPressFun);
        PTrend(i) = iP;
        ithermExp = igam*iCv./(iKT.*Vcurr);
        thermExpTrend(i) = ithermExp;
        Vnext = Vcurr*(1+ithermExp*dT);
    end
end
