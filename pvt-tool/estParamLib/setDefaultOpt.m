% setDefaultOpt - set remaining default options
function opt = setDefaultOpt(opt,optDefault)
    if(isempty(opt))
        opt = optDefault;
        return;
    end
    optNameList = fieldnames(optDefault);
    for(i=1:length(optNameList))
        ioptName = optNameList{i};
        if(~isfield(opt,ioptName))
            setfield(opt,ioptName,getfield(optDefault,ioptName));
        end
    end
end
