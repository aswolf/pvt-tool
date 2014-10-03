function [pfit,applyPriorFlag] = fitParams(pinit,nLogPFun,varargin)
    minFunAlgoDefault = @fminunc;
    NFitIterDefault = 2;
    if(isempty(varargin))
        minFunAlgo = minFunAlgoDefault;
        NFitIter   = NFitIterDefault;
    elseif(length(varargin)==1)
        minFunAlgo = varargin{1};
        NFitIter   = NFitIterDefault;
    else
        minFunAlgo = varargin{1};
        NFitIter   = varargin{2};
    end

    pfit = pinit;
    for(iter=1:fitIterNum)
        [pfit] = fminunc(nLogPFun,pfit);
    end
end
