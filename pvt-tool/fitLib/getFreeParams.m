% getFreeParams - get param. subset that has nonzero priors and are not fixed
function [pinitFree,priorFree,priorcovFree] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag)
    if(isempty(fixFlag))
        fixFlag = ones(size(pinit));
    end
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
