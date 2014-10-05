% getAllParams - get full param. vec. with free params replaced by current val
function [p] = getAllParams(pFree,p,fixFlag)
    assert(length(fixFlag) == length(p),...
        'param and fixFlag arrays must have equal length.');
    assert(sum(~fixFlag) == length(pFree),...
        'Free param array length must match number of unfixed params.');

    if(isempty(fixFlag))
        p = pFree;
    else
        p(~fixFlag) = pFree;
    end
end
