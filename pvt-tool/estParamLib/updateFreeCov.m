% updateFreeCov - update cov matrix with resutls from fitting free param subset
%
% This function is used to fill in an updated covariance 
% matrix based on fitting data that constrains only a subset of the params.
%
% Assumption is that all fixed params have their prior values
function [pcov] = updateFreeCov(pcovFree,priorcov,fixFlag)
    if(isempty(fixFlag))
        fixFlag = zeros([1 size(priorcov,1)]);
    end
    checkInput(pcovFree,priorcov,fixFlag);
    if(all(~fixFlag))
        pcov = pcovFree;
        return;
    end

    indFree = find(~fixFlag);
    pcov = priorcov;
    for(i=1:length(indFree))
        for(j=1:length(indFree))
            pcov(indFree(i),indFree(j)) = pcovFree(i,j);
        end
    end
end
function checkInput(pcovFree,priorcov,fixFlag)
    assert(all(size(priorcov)==length(fixFlag)*[1 1]),...
        'priorcov must be a square matrix with dimensions equal to fixFlag len.');
    assert(all(size(pcovFree)==sum(~fixFlag)*[1 1]),...
        ['pcovFree must be a square matrix with dimensions equal to all '...
        'unfixed params.']);
end
