function pEosDraw = drawRandEos(eosMod,Ndraw)
    if(isempty(Ndraw))
        Ndraw = 1;
    end
    pEos    = eosMod.pEos;
    pEosCov = eosMod.pEosCov;
    mask = ~isnan(diag(pEosCov)') & ~isinf(diag(pEosCov)');
    pEosCovMask = pEosCov(mask,mask);
    if(isempty(pEosCov))
        pEosDraw = repmat(pEos,Ndraw,1);
    else
        pEosDrawMask = mvnrnd(pEos(mask),pEosCovMask,Ndraw);
        pEosDraw = zeros(Ndraw,length(pEos));
        pEosDraw(:,mask)  = pEosDrawMask;
        pEosDraw(:,~mask) = repmat(pEos(~mask),Ndraw,1);
    end
end
