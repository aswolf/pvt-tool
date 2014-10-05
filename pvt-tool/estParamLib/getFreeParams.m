% getFreeParams - get param. subset that has nonzero priors and are not fixed
function [pinitFree,priorFree,priorcovFree] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag)
    if(isempty(fixFlag))
        fixFlag = ones(size(pinit));
    end

    checkInput(pinit,prior,priorcov,fixFlag);
    
    isFixed = fixFlag(:) | diag(priorcov)==0;
    indFix  = find(isFixed);
    indFree = find(~isFixed);

    pinitFree = pinit(indFree);
    priorFree = prior(indFree);
    priorcovFree = zeros(length(indFree));

    for(i=1:length(indFree))
        for(j=1:length(indFree))
            priorcovFree(i,j) = priorcov(indFree(i),indFree(j));
        end
    end
end
function checkInput(pinit,prior,priorcov,fixFlag)
    assert(length(pinit) == length(prior),...
        'pinit and prior arrays must have equal length.');
    assert(length(fixFlag) == length(pinit),...
        'pinit and fixFlag arrays must have equal length.');

    assert(all(size(priorcov) == length(prior)*[1 1]),...
        'priorcov must be square matrix with dim. equal to prior length.');

    assert(sum(~fixFlag) == length(pinit),...
        'Free param array length must match number of unfixed params.');
end
