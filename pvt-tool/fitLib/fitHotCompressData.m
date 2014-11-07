% fitHotCompressData - fit hot compression data with Eos
function [pfit pfitcov nLogPFun PressTotFun opt] = fitHotCompressData(pinitEos,fixFlag,...
        T0,NpCold,priorEos,priorcovEos,coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun,P,V,T,PerrTot,opt)

    TOLT0 = .03;
    assert(~all(fixFlag),'Nothing to fit edge case. Need to fix.')
    % Determine which measurements are 'cold' (at reference temp)
    indT0 = find(abs((T-T0)/T0) < TOLT0);
    indT  = find(abs((T-T0)/T0) >= TOLT0);

    % May need to add eos fitting specific options
    %   getFitEosDefaultOpt
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    if(isempty(fixFlag))
        fixFlag = 0;
    end
    if(length(fixFlag)==1)
        fixFlag = fixFlag*ones(size(pinitEos));
    end
    if(length(PerrTot)==1)
        PerrTot = PerrTot*ones(size(V));
    end

    %fixFlagColdOnly = fixFlag;
    %fixFlagColdOnly(NpCold+1:end) = 1;
    %[pinitCold,priorCold,priorcovCold] = ...
    %    getFreeParams(pinitEos,priorEos,priorcovEos,fixFlagColdOnly);
    pinitCold = pinitEos(1:NpCold);
    priorCold = priorEos(1:NpCold);
    priorcovCold = priorcovEos(1:NpCold,1:NpCold);
    fixFlagCold = fixFlag(1:NpCold);


    % Fit Cold subset of data
    [pfitCold pfitcovCold] = fitColdCompressData(pinitCold,fixFlagCold,...
        priorCold,priorcovCold,coldEosFun,P(indT0),V(indT0),PerrTot(indT0),opt);

    %update prior with cold fit
    fixHotEosFlag = zeros(size(fixFlag));
    fixHotEosFlag(NpCold+1:end) = 1;
    [pfit] = getAllParams(pfitCold,pinitEos,fixHotEosFlag);
    [pfitcov] = updateFreeCov(pfitcovCold,priorcovEos,fixHotEosFlag);

    PressTotFun = @(V,T,pEos)(calcPressThermAddEos(V,T,T0,...
        pEos(1:NpCold),pEos(NpCold+1:end),coldEosFun,hotEosFun,hotExtraInputs,...
        addedThermPressFun));

    wtResidFun = @(pEos)( (P(indT)-PressTotFun(V(indT),T(indT),pEos)) ...
        ./ PerrTot(indT) );

    %fixFlagHot0 = fixFlag;
    %fixFlagHot0(1:NpCold) = 1;
    %pfitHot0 = pfit;
    %[pfit nLogPFun] = fitFreeParamsWithPrior(pfit,fixFlagHot0,...
    %    pfit,pfitcov,wtResidFun,opt);

    [pfit nLogPFun] = fitFreeParamsWithPrior(pfit,fixFlag,...
        pfit,pfitcov,wtResidFun,opt);

    [pfitcov] = estParamCov(nLogPFun,pfit,fixFlag,opt);
end
%
%    pinitColdEos = pinitEos(1:NpCold);
%    pinitHotEos  = pinitEos(NpCold+1:end);
%
%    priorColdEos = priorEos(1:NpCold);
%    priorHotEos  = priorEos(NpCold+1:end);
%
%    priorerrEos = sqrt(diag(priorcovEos))';
%    priorcorrEos = priorcovEos./(priorerrEos(:)*priorerrEos(:)');
%    indInfVar = find(isinf(priorerrEos)|priorerrEos==0);
%    priorcorrEos(indInfVar,indInfVar) = 0;
%    priorcorrEos(sub2ind(size(priorcovEos),indInfVar,indInfVar)) = 1;
%    priorcovColdEos = priorcorrEos(1:NpCold,1:NpCold)...
%        .*(priorerrEos(1:NpCold)*priorerrEos(1:NpCold)');
%
%    pfiterrCold = sqrt(diag(pfitcovCold));
%    pfitcorrCold = pfitcovCold./(pfiterrCold(:)*pfiterrCold(:)');
%
%    % Update initial point and prior according to cold fit
%    pinitEos(1:NpCold)              = pfitCold;
%    priorEos(1:NpCold)              = pfitCold;
%    priorcorrEos(1:NpCold,1:NpCold) = pfitcorrCold;
%    priorerrEos(1:NpCold)           = pfiterrCold;
%
%    priorcovEos = priorcorrEos.*(priorerrEos(:)*priorerrEos(:)');
%
