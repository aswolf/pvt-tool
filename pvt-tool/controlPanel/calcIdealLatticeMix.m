function eosModIdeal = calcIdealLatticeMix(Xideal, Xendmem, eosModEndmem,...
        Pbnds,Tbnds)
    assert(length(Xendmem)==length(eosModEndmem),...
        'Must provide comp for each eosMod.');
    assert(length(Xendmem)==2, 'Only defined for 2 endmembers')
    Nendmem = 2;
    Ngrid = 10;
    T0 = 300;
    NpCold = 3;
    coldEosFun = eosModEndmem(1).coldEosFun;
    hotEosFun = eosModEndmem(1).hotEosFun;
    hotExtraInputs = eosModEndmem(1).hotExtraInputs;
    addedThermPressFun = eosModEndmem(1).addedThermPressFun;


    pEos1 = eosModEndmem(1).pEos;
    pEosEndmem = zeros(2,length(pEos1));
    pEosCovEndmem = zeros(2,length(pEos1),length(pEos1));
    fixFlag = zeros(size(pEos1));
    for(i=1:Nendmem)
        ipEos = eosModEndmem(i).pEos;
        ipEosCov = eosModEndmem(i).pEosCov;
        fixFlag(isnan(diag(ipEosCov))) = 1;
        pEosEndmem(i,:) = ipEos;
        pEosCovEndmem(i,:,:) = ipEosCov;
    end

    Vedges = zeros(Nendmem,2,2);
    for(i=1:Nendmem)
        ieosMod = eosModEndmem(i);
        iV0 = ipEos(1);
        evalP = @(V,T)(evalPressEos([],ieosMod,V,T));
        Vedges(i,1,1) = fzero(@(V)(evalP(V,Tbnds(1))-Pbnds(1)), iV0);
        Vedges(i,1,2) = fzero(@(V)(evalP(V,Tbnds(2))-Pbnds(1)), iV0);
        Vedges(i,2,1) = fzero(@(V)(evalP(V,Tbnds(1))-Pbnds(2)), iV0);
        Vedges(i,2,2) = fzero(@(V)(evalP(V,Tbnds(2))-Pbnds(2)), iV0);
    end
    Vbnds = [min(Vedges(:)) max(Vedges(:))];

    Vgridv = linspace(Vbnds(1),Vbnds(2),Ngrid);
    Tgridv = linspace(Tbnds(1),Tbnds(2),Ngrid);
    [Vgrid,Tgrid] = meshgrid(Vgridv,Tgridv);
    PerrTot = 0.1*ones(size(Vgrid));

    PgridEndmem = zeros(Nendmem,Ngrid,Ngrid);
    for(i=1:Nendmem)
        ieosMod = eosModEndmem(i);
        ipEos = ieosMod.pEos;
        %evalP = @(V,T)(evalPressEos(ipEos,ieosMod,V,T));
        iPgrid = evalPressEos(ipEos,ieosMod,Vgrid(:),Tgrid(:));
        PgridEndmem(i,:,:) = reshape(iPgrid(:),Ngrid,[]);
    end
    %imagesc(Vgridv,Tgridv,squeeze(PgridEndmem(1,:,:)));colorbar;
    %imagesc(Vgridv,Tgridv,squeeze(PgridEndmem(2,:,:)));colorbar;

    % Eval ideal lattice pressure (arithmetic mean at const vol)
    wtendmem = (Xideal-Xendmem(1))./(Xendmem(2)-Xendmem(1));
    PgridIdeal = zeros(length(Xideal),Ngrid,Ngrid);
    for(i=1:length(wtendmem))
        iwtendmem = wtendmem(i);
        PgridIdeal(i,:,:) = (1-iwtendmem)*squeeze(PgridEndmem(1,:,:)) + ...
            iwtendmem*squeeze(PgridEndmem(2,:,:));
        %imagesc(Vgridv,Tgridv,squeeze(PgridIdeal(i,:,:)));colorbar;
        %pause
    end

    % Loop over ideal samples and fit each with an Eos
    for(i=1:length(wtendmem))
        iwtendmem = wtendmem(i);
        iwt0 = iwtendmem;
        iwt0 = min(iwt0,1);
        iwt0 = max(iwt0,0);
        ipEos0 = (1-iwt0)*pEosEndmem(1,:)+iwt0*pEosEndmem(2,:);
        ipEosErr0 = Inf*ones(size(ipEos0));
        ipEosErr0(fixFlag==1) = NaN;
        ipEosCov0 = diag(ipEosErr0);
        iPgridIdeal = PgridIdeal(i,:,:);

        opt.fitcov = false;

        [ipfitEos ipfitcovEos inLogPFun iPressTotFun iopt] = ...
            fitHotCompressData(ipEos0,fixFlag,...
            T0,NpCold,ipEos0,ipEosCov0,coldEosFun,hotEosFun,hotExtraInputs,...
            addedThermPressFun,iPgridIdeal(:),Vgrid(:),Tgrid(:),PerrTot,opt);
        % iPressTotFun(Vgrid(:),Tgrid(:),pfitEos)-iPgridIdeal(:)
        ieosMod = eosModEndmem(1);
        ieosMod.name = 'ideal mix';
        ieosMod.material = 'MgFePv';
        ieosMod.pEos = ipfitEos;
        ieosMod.pEosCov = [];

        eosModIdeal(i) = ieosMod;
    end

    %pEosIdealFit = reshape([eosModIdeal.pEos],[],length(eosModIdeal))';
    %plot(Xideal,pEosIdealFit(:,1),'k-')
end
