% calcPressThermAddEos - calc total press using additive thermal contribution
function [P,KT,Cv,gam] = calcPressThermAddEos(V,T,T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,elecThermPressFun)

    V0 = pColdEos(1);

    [PCold,KTCold] = coldEosFun(V,pColdEos);
    [PHot,KTHot,CvHot,gamHot] = hotEosFun(V,T,V0,T0,pHotEos);

    if(isempty(elecThermPressFun))
        Pelec = 0;
        Cvelec = 0;
        dPdT_Velec = 0;
    else
        % Assumes that electronic thermal pressure is independent of volume
        [Pelec Cvelec dPdT_Velec] = elecThermPressFun(T);
        % Check whether electronic heat capacity and thermal press derivs are
        % calculated
        if(isempty(Cvelec))
            Cvelec = zeros(size(Pelec));
        end
        if(isempty(dPdT_Velec))
            dPdT_Velec = zeros(size(Pelec));
        end
    end

    dPdT_VHot = CvHot.*gamHot./V;

    % Additive quantities are derivatives of energy terms
    P  = PCold  + PHot + Pelec;
    KT = KTCold + KTHot;
    Cv = CvHot  + Cvelec;
    % NOTE: dPdT_V adds linearly, so we must combine terms and then determine
    % gamma after
    dPdT_V = dPdT_VHot + dPdT_Velec;
    gam = dPdT_V.*V./Cv;
end
