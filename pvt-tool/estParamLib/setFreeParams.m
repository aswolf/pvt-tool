% setFreeParams - get free param. vec. with free params replaced by current val
function [p,pcov] = setFreeParams(pFree,pcovFree,pinit,prior,priorcov,fixFlag)
    checkInput(pFree,pcovFree,pinit,prior,priorcov,fixFlag);

    indFree = find(~fixFlag);
    p(~fixFlag) = pFree;
    pcov = zeros(size(priorcov));
    for(i=1:length(indFree))
        %p(indFree(i)) = pFree(i);
        for(j=1:length(indFree))
            pcov(indFree(i),indFree(j)) = pcovFree(i,j);
        end
    end
end
function checkInput(pFree,pcovFree,pinit,prior,priorcov,fixFlag)
    assert(length(pinit) == length(prior),...
        'pinit and prior arrays must have equal length.');
    assert(length(fixFlag) == length(pinit),...
        'pinit and fixFlag arrays must have equal length.');

    assert(all(size(priorcov) == length(prior)*[1 1]),...
        'priorcov must be square matrix with dim. equal to prior length.');
    assert(all(size(pcovFree) == length(pFree)*[1 1]),...
        'pcovFree must be square matrix with dim. equal to pFree length.');

    assert(sum(~fixFlag) == length(pFree),...
        'Free param array length must match number of unfixed params.');
end
