% fitColdCompressData - fit isothermal compression data with Eos
function [pfit pfitcov nLogPFun] = fitColdCompressData(pInitEos,fixFlag,...
        priorEos,priorcovEos,coldEosFun,P,V,PerrTot,opt)
    % May need to add eos fitting specific options
    %   getFitEosDefaultOpt
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);
    if(isempty(fixFlag))
        fixFlag = 0;
    end
    if(length(fixFlag)==1)
        fixFlag = fixFlag*ones(size(pInitEos));
    end
    if(length(PerrTot)==1)
        PerrTot = PerrTot*ones(size(V));
    end
    
    %wtResidFun = @(pEosFree,pEosAll,fixFlag)( (P-...
    %    coldEosFun(V,getAllParams(pEosFree,pEosAll,fixFlag))) ...
    %    ./ PerrTot );
    wtResidFun = @(pEos)( (P-coldEosFun(V,pEos)) ...
        ./ PerrTot );

    [pfit nLogPFun] = fitFreeParamsWithPrior(pInitEos,fixFlag,...
        priorEos,priorcovEos,wtResidFun,opt);

    [pfitcov] = estParamCov(nLogPFun,pfit,fixFlag,opt);
end
