function [pcov] = corrToCov(perr,pcorr)
    pcorrT = pcorr';
    assert(all(pcorr(:) == pcorrT(:)),'corr matrix is invalid because it is not symmetric');
    eigcorr = eig(pcorr);
    assert( all(eigcorr>= 0),'corr matrix is invalid because it is not symmetric');

    perrM = perr(:)*perr(:)';
    perrM(isinf(perrM)) = 0;

    % Ensure that correlated variables with unconstrained ones don't become
    % unconstrained
    infInd = find(isinf(perr));
    perrM(sub2ind(size(perrM),infInd,infInd)) = Inf;
    pcov = perrM.*pcorr;

    %% Ensure no NaNs due to 0*Inf
    pcov(isnan(pcov)) = 0;

    %% Ensure Inf covar if either row or col is Inf
    %for(i=1:length(perr))
    %    for(j=1:length(perr))
    %        if(isinf(perr(i))|isinf(perr(j)))
    %            perrM(i,j) = Inf;
    %        end
    %    end
    %end

    %% Scale cov matrix by error matrix
    %pcorr = pcov./perrM;


    %% Ensure diagonals are 1 (by definition)
    %pcorr(sub2ind(size(pcorr),[1:length(perr)],[1:length(perr)])) = 1;
end
