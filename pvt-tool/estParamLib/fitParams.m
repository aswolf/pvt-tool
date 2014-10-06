% fitParams - fit parameters using neg. logP function w/ additional options
function [pfit] = fitParams(pinit,nLogPFun,opt)
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    pfit = pinit;
    NfitIter  = opt.NfitIter;
    minAlgFun = opt.minAlgFun;
    for(iter=1:NfitIter)
        [pfit] = minAlgFun(nLogPFun,pfit,opt);
    end
end
