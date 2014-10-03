function nLogPFun = makeNLogPFun(wtResidFun,prior,priorcov,...
        doRobustFit,robustNormParam)

    if(isempty(robustNormParam))
        robustNormParam = 5;
    end

    % NOTE: robust refers only data fitting, not prior
    %   - maybe it should apply to both?


    % Filter out unconstrained variables from prior weighting
    %  - unconstrained variables shown by infinite variance
    infVar = isinf(diag(priorcov));
    priorcovSub = priorcov(~infVar,~infVar);
    priorHessSub = inv(priorcovSub);
    priorSub = prior(~infVar);
    applyPriorFlag = ~infVar;

    priorWt = @(ptest)(0.5*(ptest(~infVar)-priorSub)...
        *priorHessSub*(ptest(~infVar)-priorSub)');

    % Construct total posterior Function nLogPFun
    %   -depends on whether it is robust or not
    if(doRobustFit)
        nu = robustNormParam;
        nLogPFun = @(ptest)(0.5*(nu+1)*sum(log(1+wtResidFun(ptest).^2/nu)) + ...
            priorWt(ptest));
    else
        nLogPFun = @(ptest)(0.5*sum(wtResidFun(ptest).^2) + ...
            priorWt(ptest));
    end
end
