function V = invertPressEos(P,T,pColdEos,pHotEos,T0,coldEosFun,hotEosFun,...
        hotExtraInputs,addedThermPressFun)

    if(numel(T)==1)
        T1 = T(1);
        T = T1*ones(numel(T),1);
    end

    V = zeros(size(P));
    V0 = pColdEos(1);
    for(i=1:length(P))
        iP = P(i);
        iT = T(i);
        PresidFun = @(ilogVf)(iP-calcPressThermAddEos(exp(ilogVf)*V0,iT,T0,...
            pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
            addedThermPressFun));

        ilogVf0 = -0.05;
        ilogVf  = fzero(PresidFun,ilogVf0);

        iV = exp(ilogVf)*V0;

        V(i) = iV;
    end
end
