% getFreeParams - get param. subset that has nonzero priors and are not fixed
function [pinitFree,priorFree,priorcovFree] = ...
        getFreeParams(pinit,prior,priorcov,fixFlag)
    if(isempty(fixFlag))
        fixFlag = zeros(size(pinit));
    end

    checkInput(pinit,prior,priorcov,fixFlag);
    
    %update fixFlag to include zero values on priorcov diag
    %fixFlag = fixFlag(:) | diag(priorcov)==0;
    indFix  = find(fixFlag);
    indFree = find(~fixFlag);

    pinitFree = pinit(indFree);
    priorFree = prior(indFree);
    priorcovFree = zeros(length(indFree));

    for(i=1:length(indFree))
        for(j=1:length(indFree))
            priorcovFree(i,j) = priorcov(indFree(i),indFree(j));
        end
    end

    checkOutput(pinitFree,priorFree,priorcovFree,fixFlag);
end
function checkInput(pinit,prior,priorcov,fixFlag)
    assert(length(pinit) == length(prior),...
        'pinit and prior arrays must have equal length.');
    assert(length(fixFlag) == length(pinit),...
        'pinit and fixFlag arrays must have equal length.');

    assert(all(size(priorcov) == length(prior)*[1 1]),...
        'priorcov must be square matrix with dim. equal to prior length.');
    diagFree = sub2ind(size(priorcov),find(fixFlag==0),find(fixFlag==0));
    assert(all(priorcov(diagFree)~=0),...
        ['Diagonal elements of cov matrix cannot be zero unless '...
        'they are fixed params. Use fixFlag to fix params']);

end
function checkOutput(pinitFree,priorFree,priorcovFree,fixFlag)
    assert(sum(~fixFlag) == length(pinitFree),...
        'Free param array length must match number of unfixed params.');
    assert(length(priorFree) == length(pinitFree),...
        'Free param arrays for init and prior must match length.');
    assert(all(size(priorcovFree) == length(priorFree)*[1 1]),...
        'priorcovFree must be square matrix with dim. equal to priorFree length.');
    assert(all(diag(priorcovFree)~=0),...
        ['Diagonal elements of cov matrix cannot be zero. '... 
        'Should not be possible for free parameters. Needs immediate debugging.']);
end
