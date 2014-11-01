function V = invertPressEos(P,T,pColdEos,pHotEos,T0,coldEosFun,hotEosFun,...
        hotExtraInputs,addedThermPressFun)

    RELTOLV = 1e-6;
    if(numel(T)==1)
        T1 = T(1);
        T = T1*ones(numel(T),1);
    end

    V = zeros(size(P));
    V0 = pColdEos(1);
    ilogVf0 =  -.05;

    for(i=1:length(P))
        iP = P(i);
        iT = T(i);

        idellogVf = Inf;
        ilogVfnext = ilogVf0;
        while(abs(idellogVf) > RELTOLV)
            ilogVf = ilogVfnext;
            [iPcalc,iKT] = calcPressThermAddEos(exp(ilogVf)*V0,iT,T0,...
                pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
                addedThermPressFun);
            idellogVf = (iPcalc-iP)/iKT;
            ilogVfnext = ilogVf + idellogVf;
        end


        iV = exp(ilogVf)*V0;

        V(i) = iV;
    end
end
