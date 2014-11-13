% fitErrModResid - fit error bar adjustment params given current model residuals
function [pfitErrModList nLogPFunList] = fitErrModResidMultiGrp(measGrpID,...
        pinitErrModList,priorErrModList,priorcovErrModList,linTransM,...
        yresid,dydxmod,xerr,opt)

    % May need to add eos fitting specific options
    optDefault = getEstParamDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    % If unspecified, measGrpID = '1' for all datapoints
    if(isempty(measGrpID))
        measGrpID = '1';
    end
    % If only one group specified, measGrpID is the same for all datapoints
    if(isstr(measGrpID))
        Ndat = length(yresid);
        measGrpID = cellstr(repmat(measGrpID,Ndat,1));
    end
    uniqID = unique(measGrpID);
    NmeasGrp = length(uniqID);

    if(NmeasGrp==1)
        [pfitErrModList nLogPFunList] = fitErrModResid(pinitErrModList,...
            priorErrModList,priorcovErrModList,linTransM, ...
            yresid,dydxmod,xerr,opt);
        return
    end



    checkInput(measGrpID,pinitErrModList,priorErrModList,priorcovErrModList,...
        linTransM,yresid,dydxmod,xerr);

    pfitErrModList = zeros(size(pinitErrModList));

    for(i=1:NmeasGrp)
        ipinitErrMod    = pinitErrModList(i,:);
        ipriorErrMod    = priorErrModList(i,:);
        ipriorcovErrMod = squeeze(priorcovErrModList(i,:,:));
        imeasGrpInd = find(strcmp(measGrpID,uniqID{i}));
        iyresid = yresid(imeasGrpInd);
        idydxmod = dydxmod(imeasGrpInd,:);
        ixerr = xerr(imeasGrpInd,:);
        
        [ipfitErrMod inLogPFun] = fitErrModResid(ipinitErrMod,...
            ipriorErrMod,ipriorcovErrMod,linTransM,iyresid,idydxmod,ixerr,opt);
        pfitErrModList(i,:) = ipfitErrMod;
        nLogPFunList{i}     = inLogPFun;
    end
end

function checkInput(measGrpID,pinitErrModList,priorErrModList,...
        priorcovErrModList,linTransM,yresid,dydxmod,xerr)
    Ndat = size(yresid,1);
    dimx = size(dydxmod,2);
    uniqID = unique(measGrpID);
    NmeasGrp = length(uniqID);
    NpErrMod = size(pinitErrModList,2);

    assert(ndims(pinitErrModList)==2,'pinitErrModList must be 2 dimensional')
    assert(ndims(priorErrModList)==2,'priorErrModList must be 2 dimensional')

    covErrModGrp = squeeze(priorcovErrModList(1,:,:));
    assert(size(covErrModGrp,1)==size(covErrModGrp,2),...
        'priorcovErrModList must be a list of square matrices')
    if(~isempty(linTransM))
        assert(size(linTransM,1)==dimx,...
            'linTransM must have one row per x dimension.');
        assert(size(linTransM,2)==NpErrMod+1,['linTransM must have one col per '...
            'errMod parameter plus one extra for offsets.']);
    end

    assert(size(pinitErrModList,1)==NmeasGrp,...
        'pinitErrModList must have len NmeasGrp')
    assert(size(priorErrModList,1)==NmeasGrp,...
        'priorErrModList must have len NmeasGrp')

    assert(iscellstr(measGrpID),...
        'measGrpID must be a cell array of strings identifying each group');
    assert(length(measGrpID)==Ndat,'measGrpID must have one entry per datum.')
    assert(size(pinitErrModList,1)==NmeasGrp,...
        'The errModel must have a set of parameters for each measGrp'); 

    %NpErrMod = size(pinitErrModList,2);
    %assert(all(size(priorErrMod)==size(pinitErrMod)),...
    %    'Number of error mod params must be equal in prior and init');
    %assert(all(size(priorcovErrModList)==[NmeasGrp NpErrMod*[1 1]]),...
    %    'priorcov matrix must provide square with dimensions of param num');
end
