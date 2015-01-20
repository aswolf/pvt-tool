function [Pad,Vad,Tad,Kad,gamad,thmExpad]=calcAdiabat(Tfoot,Pfoot,Pstop,dlogV,...
        eosMod,pEos)
    % pEos=[] -> use best-fit
    if(~exist('pEos'))
        pEos = [];
    end
    if(isempty(pEos))
        pEos = eosMod.pEos;
    end


    evalP = @(V,T)(evalPressEos(pEos,eosMod,V,T));
    V0 = eosMod.pEos(1);
    Vfoot = fzero(@(V)(evalP(V,Tfoot)-Pfoot),V0);

    dPfoot = Pstop-Pfoot;
    if(dPfoot < 0)
        dlogV = abs(dlogV);
    end

    iV = Vfoot;
    iT = Tfoot;
    ind=1;
    while(true)
        [iP,iPderivs,iKT,iCv,igam,ithmExp] = evalP(iV, iT);
        iKad = iKT*(1+ithmExp*igam*iT);
        Pad(ind) = iP;
        Vad(ind) = iV;
        Tad(ind) = iT;
        Kad(ind) = iKad;
        gamad(ind)= igam;
        thmExpad(ind) = ithmExp;

        if(abs(iP-Pfoot) > abs(dPfoot))
            break;
        else
            iT = iT*(1-igam*dlogV);
            iV = iV*(1+dlogV);
            ind = ind+1;
        end
    end
end
