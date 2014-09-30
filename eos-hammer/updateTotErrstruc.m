function datStruc = updateTotErrstruc(eosStruc,datStruc,doRobustFit,updatePscl)
    clf;
    hold on;
    for(i=1:length(datStruc))
        idatStruc = updateTotErrstruc1(eosStruc,datStruc(i),doRobustFit,updatePscl);
        datStruc(i) = idatStruc;
    end
    hold off;
end
function datStruc = updateTotErrstruc1(eosStruc,datStruc,doRobustFit,updatePscl)
    VfracerrDefault = 3e-3;
    % Determine scale factor for relative Vol errors using residual scatter of
    % data
    runID    = datStruc.runID    ;
    P        = datStruc.P        ;
    Perr     = datStruc.Perr     ;
    T        = datStruc.T        ;
    Terr     = datStruc.Terr     ;
    V        = datStruc.V        ;
    Verr     = datStruc.Verr     ;
    Vmark    = datStruc.Vmark    ;
    Vmarkerr = datStruc.Vmarkerr ;
    PerrTot  = datStruc.PerrTot  ;

    if(any(Vmarkerr==0))
        Vmarkerr = VfracerrDefault*Vmark;
    end


    peos = eosStruc.eos;
    eosColdFun = eosStruc.coldFun;
    eosMGDFun  = eosStruc.MGDFun;

    peosmark = eosStruc.eosmark;
    eosmarkColdFun = eosStruc.markcoldFun;
    eosmarkMGDFun  = eosStruc.markMGDFun;

    T0         = eosStruc.T0;

    PVTmod = evalEosPVTwrap(V,T,peos,T0,eosColdFun,eosMGDFun);
    PVTmark = evalEosPVTwrap(Vmark,T,peosmark,T0,eosmarkColdFun,eosmarkMGDFun);
    if(updatePscl)
        P = PVTmark.P;
        datStruc.P = P;
    end
    
    %clf;hold on;
    %scatter(P,V,50,T,'o')
    %scatter(PVTmod.P,V,50,T,'x')
    %hold off;
    scatter(P,P-PVTmod.P,50,T,'o')

    uniqID = unique(runID);

    %NEED TO do this for each uniqID!

    resid = PVTmod.P-P;
    N = length(P);
    %PvarTot0 = (PVTmod.dPdV.*Verr).^2 +  ((PVTmark.dPdT-PVTmod.dPdT).*Terr).^2;

    errModpriorwid = [.1 .1 .03];
    %logErrFacinit = zeros(1,length(uniqID)*3);
    logErrFacinitM = datStruc.logErrModFac;
    logErrFacinit = logErrFacinitM(:)';

    %nLogP = nLogPErrMod(logErrFac,runID,Verr,Vmarkerr,Terr,...
    %    PVTmod,PVTmark,doRobustFit,ErrModpriorwid);
    nLogPFun = @(logErrFac)(nLogPErrMod(logErrFac,resid,runID,Verr,Vmarkerr,Terr,...
        PVTmod,PVTmark,doRobustFit,errModpriorwid));

    varTotFun = @(logErrFac)(calcPvarTot(logErrFac,runID,Verr,Vmarkerr,Terr,PVTmod,PVTmark));

    %PvarTot = varTotFun(logErrFacinit);
    %nLogPFun(logErrFacinit)

    logErrFacFit = fminunc(nLogPFun,logErrFacinit);
    logErrFacFitM = reshape(logErrFacFit(:),[],length(uniqID))';
    PerrTot = sqrt(varTotFun(logErrFacFit));
    %ploterr([1:N],resid,[],PerrTot,'ro')
    %plot([1:N],resid./PerrTot,'ko')
    %hist(resid./PerrTot,20)

    datStruc.logErrModFac = logErrFacFitM;
    datStruc.PerrTot = PerrTot;

    %sqrt(mean((resid./sqrt(PvarTot0)).^2))
    %logVmarkerrFac = linspace(-2,2,100);
    %typResidScl = zeros(size(logVmarkerrFac));
    %for(i=1:length(logVmarkerrFac))
    %    iVmarkerrFac = exp10(logVmarkerrFac(i));
    %    iPerrTot = sqrt(PvarTot0 + (iVmarkerrFac*PVTmark.dPdV.*Vmarkerr).^2);

    %    %NOTE: Should be made robust!
    %    % Do linear regression on volume and Temp to partially account for 
    %    % systematic offset


    %    residScl = resid./iPerrTot;

    %    %designM = [T/300-1,1-V/peos(1),ones(size(V))];
    %    %designMScl = designM./repmat(iPerrTot,1,size(designM,2));
    %    %residLinMod = designM\resid;
    %    %residLinModScl = designMScl\residScl;
    %    
    %    %ploterr([1:N],resid-designM*residLinMod,[],iPerrTot,'ko')
    %    %ploterr([1:N],(resid-designM*residLinMod)./iPerrTot,[],ones(size(P)),'ko')
    %    %ploterr([1:N],(residScl-designMScl*residLinModScl),[],ones(size(P)),'ro')
    %    %typResidScl(i) = 1.4826*median(abs(residScl-designMScl*residLinModScl));
    %    typResidScl(i) = 1.4826*median(abs(residScl));
    %end
    %indLastMatch = find(log10(typResidScl)>0,1,'last');

    %VmarkerrFac = exp10(interp1(log10(typResidScl(indLastMatch+[0 1])),...
    %    logVmarkerrFac(indLastMatch+[0 1]),0));

    %
    %PerrTot = sqrt(PvarTot0 + (VmarkerrFac*PVTmark.dPdV.*Vmarkerr).^2);

end
function PvarTot = calcPvarTot(logErrFac,runID,Verr,Vmarkerr,Terr,PVTmod,PVTmark)
    uniqID = unique(runID);
    NID = length(uniqID);
    logErrFacM = reshape(logErrFac(:),[],NID)';
    PvarTot = zeros(size(Verr));
    
    for(i=1:length(uniqID))
        iInd = find(runID==uniqID(i));
        ilogErrFac = logErrFacM(i,:);

        iPvarTot = (exp10(ilogErrFac(1))*PVTmod.dPdV(iInd).*Verr(iInd)).^2+...
            (exp10(ilogErrFac(2))*PVTmark.dPdV(iInd).*Vmarkerr(iInd)).^2+...
            (exp10(ilogErrFac(3))*(PVTmark.dPdT(iInd)-PVTmod.dPdT(iInd)).*Terr(iInd)).^2;
        PvarTot(iInd) = iPvarTot;
    end
end

function nLogP = nLogPErrMod(logErrFac,resid,runID,Verr,Vmarkerr,Terr,...
        PVTmod,PVTmark,varargin)

    errModpriorwidDefault = [.1 .1 .03];
    if(isempty(varargin))
        doRobustFit = false;
        errModpriorwid = errModpriorwidDefault;
    elseif(length(varargin)==1)
        doRobustFit = varargin{1};
        errModpriorwid = errModpriorwidDefault;
    else
        doRobustFit = varargin{1};
        errModpriorwid = varargin{2};
    end

    VDiffpriorwid   = errModpriorwid(1);
    Tpriorwid       = errModpriorwid(2);
    runDiffpriorwid = errModpriorwid(3);

    uniqID = unique(runID);
    NID = length(uniqID);
    logErrFacM = reshape(logErrFac(:),[],NID)';

    PvarTot = calcPvarTot(logErrFac,runID,Verr,Vmarkerr,Terr,...
        PVTmod,PVTmark);
    if(doRobustFit)
        nu = 5;
        nLogPErrMod=0.5*sum((nu+1)*log(1+resid.^2./PvarTot/nu)+log(PvarTot));
    else
        nLogPErrMod=0.5*sum(resid.^2./PvarTot+log(PvarTot));
    end

    nLogPrior = 0;
    IDnum = [1:length(uniqID)];
    for(i=1:length(uniqID))
        iIDother = IDnum~=i;
        ilogErrFac = logErrFacM(i,:);

        irunDiff = logErrFacM(iIDother,:)-ones(sum(iIDother),1)*ilogErrFac;
        irunDiffCost = 0.5*sum((irunDiff(:)./runDiffpriorwid).^2);
        iVDiffCost = 0.5*(diff(ilogErrFac(1:2))/VDiffpriorwid)^2;
        iTCost = 0.5*(ilogErrFac(3)/Tpriorwid)^2;
        nLogPrior = nLogPrior + irunDiffCost + iVDiffCost + iTCost;
    end
    nLogP = nLogPErrMod + nLogPrior;
end

