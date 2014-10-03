% getAllParams - get full param. vec. with free params replaced by current val
function [p] = getAllParams(pFree,p,fixFlag)
    if(isempty(fixFlag))
        p = pFree;
    else
        p(~fixFlag) = pFree;
    end
end
