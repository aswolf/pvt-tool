% calcPressThermAddEos - calc total press using additive thermal contribution
function [P,KT,Cv,gam,thmExp] = calcPressThermAddEos(V,T,T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,hotExtraInputs,elecThermPressFun)
    V0 = pColdEos(1);

    [PCold,KTCold] = coldEosFun(V,pColdEos);
    [PHot,KTHot,CvHot,gamHot] = hotEosFun(V,T,V0,T0,pHotEos,hotExtraInputs{:});

    if(isempty(elecThermPressFun))
        Pelec = 0;
        Cvelec = 0;
        dPdT_Velec = 0;
    else
        % Assumes that added thermal pressure is independent of volume
        %    -> relevant for electronic pressure contribution in metals
        [Pelec Cvelec dPdT_Velec] = elecThermPressFun(T);
        % Check whether added heat capacity and thermal press derivs are
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
    thmExp = gam.*Cv./(KT.*V);
end
