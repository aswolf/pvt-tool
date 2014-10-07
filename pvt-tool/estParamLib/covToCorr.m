function [perr,pcorr] = covToCorr(pcov)
    pcovT = pcov';
    assert(all(pcov(:) == pcovT(:)),'cov matrix is invalid because it is not symmetric');
    indFin = ~isinf(diag(pcov));
    eigcov = eig(pcov(indFin,indFin));
    assert( all(eigcov>= 0),'cov matrix is invalid because it is not symmetric');

    perr  = sqrt(diag(pcov))';
    perrM = perr(:)*perr(:)';

    % Ensure Inf covar if either row or col is Inf
    for(i=1:length(perr))
        for(j=1:length(perr))
            if(isinf(perr(i))|isinf(perr(j)))
                perrM(i,j) = Inf;
            end
        end
    end

    % Scale cov matrix by error matrix
    pcorr = pcov./perrM;

    % Ensure no NaNs due to 0/0
    pcorr(isnan(pcorr)) = 0;

    % Ensure diagonals are 1 (by definition)
    pcorr(sub2ind(size(pcorr),[1:length(perr)],[1:length(perr)])) = 1;
end
