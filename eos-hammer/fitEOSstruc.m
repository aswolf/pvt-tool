function [eosFitStruc jointDatStruc] = fitEOSstruc(eosInitStruc,eosPriorStruc,datStruc,...
        varargin)

    T0TolDefault = 3e-2;
    if(isempty(varargin))
        doRobustFit = false;
        T0Tol = T0TolDefault;
    elseif(length(varargin)==1)
        doRobustFit = varargin{1};
        T0Tol = T0TolDefault;
    elseif(length(varargin)==2)
        doRobustFit = varargin{1};
        T0Tol = varargin{2};
    end

    [jointDatStruc] = joinEOSDataSets(datStruc);

    P        = jointDatStruc.P        ;
    V        = jointDatStruc.V        ;
    T        = jointDatStruc.T        ;
    PerrTot  = jointDatStruc.PerrTot  ;

    eosFitStruc = eosInitStruc;

    eosColdFun = eosInitStruc.coldFun;
    eosMGDFun  = eosInitStruc.MGDFun;
    T0         = eosInitStruc.T0;

    pinit = eosInitStruc.eos;
    prior = eosPriorStruc.eos;
    priorcov = eosPriorStruc.eoscov;

    %indT0  = find(T <= T0);
    %indHot = find(T >  T0);
    indT0  = find(abs(T/T0-1) <= T0Tol);
    indHot = find(abs(T/T0-1) >  T0Tol);

    residFun = @(peosFree,peosAll,fixFlag,ind)( (P(ind)-...
        extractfield(evalEosPVTwrap(V(ind),T(ind),...
        getAllParams(peosFree,prior,fixFlag),T0,eosColdFun,eosMGDFun),'P')')...
        ./ PerrTot(ind) );


    % First fit only cold params for ambient temp measurements
    fixFlag = ones(size(prior));
    fixFlag(1:3) = 0;
    [priorFree,priorcovFree] = getFreeParams(prior,priorcov,fixFlag);
    pinitFree = pinit(fixFlag==0);

    %P(indT0)
    %getAllParams(pinitFree,prior,fixFlag)
    %modstr=evalEosPVTwrap(V(indT0),T(indT0),...
    %    getAllParams(pinitFree,prior,fixFlag),T0,eosColdFun,eosMGDFun)

    [pcoldfit pcoldcov] = fitEOSparams(pinitFree,priorFree,priorcovFree,@(pcold)residFun(pcold,prior,fixFlag,indT0),doRobustFit);
    %residFun(pcoldfit,prior,fixFlag,indT0)
    %residFun(prior,prior,zeros(size(priorHot)),indT0)

    % Translate back to full parameter set as prior for all parameter Hot data
    % fit
    [priorHot,priorcovHot] = setFreeParams(pcoldfit,pcoldcov,prior,priorcov,fixFlag);
    %residFun(priorHot,priorHot,zeros(size(priorHot)),indHot)

    fixFlagHot = (diag(priorcovHot)==0)';
    [priorHotFree,priorcovHotFree] = getFreeParams(priorHot,priorcovHot,fixFlagHot);
    %pinitHotFree = priorHot(fixFlagHot==0);
    [pfit, pfitcov] = fitEOSparams(priorHotFree,priorHotFree,priorcovHotFree,...
        @(peos)residFun(peos,priorHot,fixFlagHot,indHot),doRobustFit);

    resid = PerrTot.*residFun(pfit,priorHot,fixFlagHot,[1:length(P)]);
    jointDatStruc.residP = resid;

    % Translate back to full parameter set as prior for all parameter Hot data
    % fit
    [pfitFull,pfitcovFull] = setFreeParams(pfit,pfitcov,priorHot,priorcovHot,fixFlagHot);
    eosFitStruc.eos    = pfitFull;
    eosFitStruc.eoscov = pfitcovFull;
end
function [pfit pfitcov] = fitEOSparams(pinit,prior,priorcov,residFun,doRobustFit)
    % First filter out infinite variance (unconstrained) variables from the
    %  prior weighting
    infVar = isinf(diag(priorcov));
    priorcovSub = priorcov(~infVar,~infVar);
    priorHessSub = inv(priorcovSub);
    priorSub = prior(~infVar);
    priorWt = @(pfit)(0.5*(pfit(~infVar)-priorSub)*priorHessSub*(pfit(~infVar)-priorSub)');

    % Construct total posterior Function (nLogLkFun), according to whether it
    % is robust or not
    if(doRobustFit)
        nu = 5;
        nLogLkFun = @(pfit)(0.5*(nu+1)*sum(log(1+residFun(pfit).^2/nu)) + ...
            priorWt(pfit));
    else
        nLogLkFun = @(pfit)(0.5*sum(residFun(pfit).^2) + priorWt(pfit));
    end

    %NOTE: The hessian calculation of fminunc is NOT WORKING, so we use
    %  fastminfunc instead!
    %    -> Refit once to be sure at minimum
    [pfit,nLogLkMin] = fminunc(nLogLkFun,pinit);
    [pfitF,nLogLkMinF] = fminunc(nLogLkFun,pfit);
    [pfitcovF,scale,H,scaleFun,pfitFscl] = calcHessScale(nLogLkFun,pfitF);
    pfitF = pfitFscl;

    pfit = pfitF;
    pfitcov = pfitcovF;

    %[pfit,costval,hessOut,output] = fastminfunc(nLogLkFun,pinit);
    %[pfitF,costvalF,hessOutF,outputF] = fastminfunc(nLogLkFun,pfit,hessOut);
    %pfit = pfitF;
    %pfitcov = inv(hessOutF);
    %pfit = pfitF;
    %pfitcov = inv(hessOutF);
end
function [peosFree,peoscovFree] = getFreeParams(peos,peoscov,fixFlag)
    isFixed = diag(peoscov)==0 | fixFlag(:);
    indFix  = find(isFixed);
    indFree = find(~isFixed);

    peosFree = peos(indFree);
    peoscovFree = zeros(length(indFree));

    for(i=1:length(indFree))
        for(j=1:length(indFree))
            peoscovFree(i,j) = peoscov(indFree(i),indFree(j));
        end
    end
end
function [peos,peoscov] = setFreeParams(peosFree,peoscovFree,...
        peos,peoscov,fixFlag)
    indFree = find(~fixFlag);
    peos(~fixFlag) = peosFree;
    for(i=1:length(indFree))
        %peos(indFree(i)) = peosFree(i);
        for(j=1:length(indFree))
            peoscov(indFree(i),indFree(j)) = peoscovFree(i,j);
        end
    end
end
function [peos] = getAllParams(peosFree,peos,fixFlag)
    if(isempty(fixFlag))
        peos = peosFree;
    else
        peos(~fixFlag) = peosFree;
    end
end
