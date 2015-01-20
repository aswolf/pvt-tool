function [Texp,thmExp,Vexp,Kexp,gamexp]=calcThmExp(Pconst,Tfoot,Tstop,dT,...
        eosMod,pEos)
    % pEos=[] -> use best-fit
    if(isempty(pEos))
        pEos = eosMod.pEos;
    end

    evalP = @(V,T)(evalPressEos(pEos,eosMod,V,T));
    V0 = pEos(1);
    Vfoot = fzero(@(V)(evalP(V,Tfoot)-Pconst),V0);

    dTfoot = Tstop-Tfoot;

    iV = Vfoot;
    iT = Tfoot;
    ind=1;
    while(true)
        [iP,iPderivs,iKT,iCv,igam,ithmExp] = evalP(iV, iT);
        idlogV = ithmExp*dT;
        Vexp(ind) = iV;
        Texp(ind) = iT;
        Kexp(ind) = iKT;
        gamexp(ind)= igam;
        thmExp(ind) = ithmExp;

        if(abs(iT-Tfoot) > abs(dTfoot))
            break;
        else
            iV = iV*(1+idlogV);
            iT = iT+dT;
            ind = ind+1;
        end
    end
end
