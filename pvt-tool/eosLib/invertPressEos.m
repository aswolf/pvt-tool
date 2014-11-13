function V = invertPressEos(P,T,pColdEos,pHotEos,T0,coldEosFun,hotEosFun,...
        hotExtraInputs,addedThermPressFun)

    RELTOLV = 1e-6;
    if(numel(T)==1)
        T1 = T(1);
        T = T1*ones(numel(T),1);
    end

    V = zeros(size(P));
    V0 = pColdEos(1);
    ilogVf0 =  -.05;


    logVgrid = linspace(log(.1),log(2.3),20);
    Pgrid = zeros(size(logVgrid));
    KTgrid = zeros(size(logVgrid));
    thmExpgrid = zeros(size(logVgrid));
    for(i=1:length(logVgrid))
        iV = exp(logVgrid(i))*V0;
        [iPcalc,iKT,iCv,igam,ithmExp] = calcPressThermAddEos(iV,T0,T0,...
            pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
            addedThermPressFun);
        Pgrid(i) = iPcalc;
        KTgrid(i) = iKT;
        thmExpgrid(i) = ithmExp;
    end

    % Filter out NaN values resulting from riduculously large volumes
    realInd = find(~isnan(Pgrid));
    logVgrid = logVgrid(realInd);
    Pgrid = Pgrid(realInd);
    KTgrid = KTgrid(realInd);
    thmExpgrid = thmExpgrid(realInd);
    pp = pchipd(logVgrid,Pgrid,-KTgrid);
    logVp = linspace(log(.15),log(1.15),100);
    %plot(logVgrid,Pgrid,'kx',logVp,ppval(pp,logVp),'r-')

    logVcoldInit = zeros(size(P));
    logThmExpInit = zeros(size(P));
    for(i=1:length(P))
        iP = P(i);
        ifun = @(logV)(ppval(pp,logV)-iP);
        logVinit = interp1(Pgrid,logVgrid,iP);
        ilogVcold = fzero(ifun,logVinit);
        logVcoldInit(i) = ilogVcold;
        logThmExpInit(i) = interp1(logVgrid,log(thmExpgrid),ilogVcold,'spline','extrap');
    end

    %logVhotInit = logVcoldInit + exp(logThmExpInit).*(T-T0);



    for(i=1:length(P))
        iP = P(i);
        iT = T(i);
        %ilogVinit = ilogVf0;
        ilogVinit = logVcoldInit(i);
        %ilogVinit = logVhotInit(i);

        idellogVf = Inf;
        ilogVfnext = ilogVinit;
        while(abs(idellogVf) > RELTOLV)
            ilogVf = ilogVfnext;
            [iPcalc,iKT] = calcPressThermAddEos(exp(ilogVf)*V0,iT,T0,...
                pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
                addedThermPressFun);
            idellogVf = (iPcalc-iP)/iKT;
            ilogVfnext = ilogVf + idellogVf;
        end
        iV = exp(ilogVf)*V0;
        V(i) = iV;
    end
end

    %logVgrid = linspace(log(.2),log(5),20);
    %Pgrid = zeros(size(logVgrid));
    %KTgrid = zeros(size(logVgrid));
    %thmExpgrid = zeros(size(logVgrid));
    %for(i=1:length(logVgrid))
    %    iV = exp(logVgrid(i))*V0;
    %    [iPcalc,iKT,iCv,igam,ithmExp] = calcPressThermAddEos(iV,T0,T0,...
    %        pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
    %        addedThermPressFun);
    %    Pgrid(i) = iPcalc;
    %    KTgrid(i) = iKT;
    %    thmExpgrid(i) = ithmExp;
    %end

    %%maxT
    %Tmax = max(T);
    %PgridT = zeros(size(logVgrid));
    %KTgridT = zeros(size(logVgrid));
    %thmExpgridT = zeros(size(logVgrid));
    %for(i=1:length(logVgrid))
    %    iV = exp(logVgrid(i))*V0;
    %    [iPcalc,iKT,iCv,igam,ithmExp] = calcPressThermAddEos(iV,Tmax,T0,...
    %        pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
    %        addedThermPressFun);
    %    PgridT(i) = iPcalc;
    %    KTgridT(i) = iKT;
    %    thmExpgridT(i) = ithmExp;
    %end

    %pp = pchipd(logVgrid,Pgrid,-KTgrid);
    %logVp = linspace(log(.15),log(1.15),100);
    %plot(logVgrid,Pgrid,'kx',logVp,ppval(pp,logVp),'r-')

    %logVcoldInit = zeros(size(P));
    %logThmExpInit = zeros(size(P));
    %for(i=1:length(P))
    %    iP = P(i);
    %    ifun = @(logV)(ppval(pp,logV)-iP);
    %    logVinit = interp1(Pgrid,logVgrid,iP);
    %    ilogVcold = fzero(ifun,logVinit);
    %    logVcoldInit(i) = ilogVcold;
    %    logThmExpInit(i) = interp1(logVgrid,log(thmExpgrid),ilogVcold);
    %end

    %for(i=1:length(P))
    %    iP = P(i);
    %    iT = T(i);

    %    % First refine cold volumes
    %    idellogVf = Inf;
    %    ilogVfnext = logVcoldInit(i);
    %    while(abs(idellogVf) > RELTOLV)
    %        ilogVf = ilogVfnext;
    %        [iPcalc,iKT,iCv,igam,ithmExp] = calcPressThermAddEos(exp(ilogVf)*V0,T0,T0,...
    %            pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
    %            addedThermPressFun);
    %        idellogVf = (iPcalc-iP)/iKT;
    %        ilogVfnext = ilogVf + idellogVf;
    %    end
    %    %iVcold = exp(ilogVfnext)*V0;
    %    %
    %    %Now increase to target temperature along isobar

    %    idellogVf = Inf;
    %    while(abs(idellogVf) > RELTOLV)
    %        ilogVf = ilogVfnext;
    %        [iPcalc,iKT,iCv,igam,ithmExp] = calcPressThermAddEos(exp(ilogVf)*V0,iT,T0,...
    %            pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
    %            addedThermPressFun);
    %        idellogVf = (
    %        idellogVf = (iPcalc-iP)/iKT;
    %        ilogVfnext = ilogVf + idellogVf;
    %    end

    %    iV = exp(ilogVf)*V0;

    %    V(i) = iV;
    %end

    %plot(logVcoldInit,P,'ko')

    %logVhotInit = logVcoldInit + exp(logThmExpInit).*(T-T0);
    %PInit = zeros(size(P));

    %for(i=1:length(V))
    %    iV = exp(logVhotInit(i))*V0;
    %    iT = T(i);
    %    [iPcalc,iKT,iCv,igam,ithmExp] = calcPressThermAddEos(iV,iT,T0,...
    %        pColdEos,pHotEos,coldEosFun,hotEosFun,hotExtraInputs,...
    %        addedThermPressFun);
    %    PInit(i) = iPcalc;
    %end
