% setFreeParams - get free param. vec. with free params replaced by current val
function [p,pcov] = setFreeParams(pFree,pcovFree,pinit,prior,priorcov,fixFlag)
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
function [p] = getAllParams(pFree,p,fixFlag)
    if(isempty(fixFlag))
        p = pFree;
    else
        p(~fixFlag) = pFree;
    end
end
