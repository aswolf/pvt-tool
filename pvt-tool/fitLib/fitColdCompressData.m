% fitColdCompressData - fit isothermal compression data with Eos
function [pfit pfitcov nLogPFun] = fitColdCompressData(pInitEos,fixFlag,...
        priorEos,priorcovEos,coldEosFun,P,V,PerrTot,indDat,opt)
    % May need to add eos fitting specific options
    %   getFitEosDefaultOpt
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);
    if(isempty(fixFlag))
        fixFlag = 0;
    end
    if(isempty(indDat))
        indDat = [1:length(V)]';
    end
    if(length(fixFlag)==1)
        fixFlag = fixFlag*ones(size(pInitEos));
    end
    if(length(PerrTot)==1)
        PerrTot = PerrTot*ones(size(V));
    end
    
    %wtResidFun = @(pEosFree,pEosAll,fixFlag)( (P(indDat)-...
    %    coldEosFun(V(indDat),getAllParams(pEosFree,pEosAll,fixFlag))) ...
    %    ./ PerrTot(indDat) );
    wtResidFun = @(pEos)( (P(indDat)-coldEosFun(V(indDat),pEos)) ...
        ./ PerrTot(indDat) );

    [pfit nLogPFun] = fitFreeParamsWithPrior(pInitEos,fixFlag,...
        priorEos,priorcovEos,wtResidFun,opt);

    % what about fixFlag?
    %  Need to add that into estParamCov
    [pfitcov] = estParamCov(nLogPFun,pfit,opt);
end
