function opt = viewPVTFit(PVTeval, viewmode, opt)
    angChar = char(197);

    if(~exist('viewmode'))
        viewmode = [];
    end
    if(~exist('opt'))
        opt = [];
    end

    if(isempty(viewmode))
        viewmode = 'normal';
    end
    optDefault = getDefaultOpt();
    opt = setDefaultOpt(opt,optDefault);

    fontSz = 14;
    %
    dT = 500;
    T0 = 1000;

    eosMod = PVTeval.sampEosFit;
    PVT = [PVTeval.PVTdataList.Pmark, PVTeval.PVTdataList.V, ...
        PVTeval.PVTdataList.T];
    ndat = size(PVT,1);
    indHot = find(PVT(:,3)>310);
    indCold = setdiff([1:size(PVT,1)],indHot);

    Presid = PVT(:,1)-PVTeval.Psamp;
    PErrTot = PVTeval.PVTdataList.PErrTot;
    PbndsDat = [min(PVT(:,1)) max(PVT(:,1))];
    VbndsDat = [min(PVT(:,2)) max(PVT(:,2))];
    TbndsDat = [min(PVT(:,3)) max(PVT(:,3))];

    VbndsDat = VbndsDat + .1*diff(VbndsDat)*[-1 1];

    TmaxDat = dT*round(TbndsDat(2)/dT);
    TbndsDat(2) = TmaxDat;

    if(isempty(opt.Tbnds))
        opt.Tbnds = TbndsDat;
    end
    if(isempty(opt.Vbnds))
        opt.Vbnds = VbndsDat;
    end


    Tmod = [300 T0:dT:opt.Tbnds(2)];
    Vmod = linspace(opt.Vbnds(1),opt.Vbnds(2),100)';
    if(isempty(opt.hax))
        clf;
        opt.hax = gca();
    else
        axes(opt.hax);
        cla;
    end

    switch viewmode
        case 'hist'
            nbin = max(floor(ndat/10), 10);
            xbnd = [min(Presid./PErrTot), max(Presid./PErrTot)];
            dxbin = diff(xbnd)/nbin;
            xbnd = xbnd + 0.1*dxbin*[-1 1];
            xedges = linspace(xbnd(1),xbnd(2),nbin+1)';
            xbin = xedges(1:end-1)+dxbin/2.0;
            [nbin] = histc(Presid./PErrTot,xedges);
            [nbinHot] = histc(Presid(indHot)./PErrTot(indHot),xedges);
            [nbinCold] = histc(Presid(indCold)./PErrTot(indCold),xedges);

            hbar = bar(xedges,nbin,'histc');
            hold on;
            hbarHot = bar(xedges,nbinHot,'histc');
            hbarCold = bar(xedges,nbinCold,'histc');
            hold off;
            xbnd = xlim();
            xbnd = max(abs(xbnd))*[-1 1];
            xmod = linspace(xbnd(1),xbnd(2),100)';
            hold on;
            plot(xmod,sum(nbin)*diff(xbin(1:2))/sqrt(2*pi)*exp(-0.5*xmod.^2),...
                '-','Color',.5*[1 1 1],'LineWidth',2);
            xlim(xbnd)
            hold off;
            set(hbar,'FaceColor',.95*[1 1 1],'LineWidth',2,'EdgeColor',[0 0 0])
            set(hbarCold,'FaceColor','none','LineWidth',2,'EdgeColor',[0 0 1])
            set(hbarHot,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0])
            % xlabel('Normalized Residuals [\DeltaP /\sigma ]','FontSize',fontSz);
            % ylabel('Number of Residuals','FontSize',fontSz)
            xlabel('Normalized Residuals','FontSize',fontSz);
            ylabel('Number','FontSize',fontSz)
        case 'reduced'
            Pred = evalPressEos([],eosMod,PVT(:,2),300);
            Pmodred = evalPressEos([],eosMod,Vmod,300);

            plot(Pmodred, Vmod,'b-','LineWidth',2);
            hold on;
            scatter(PVT(:,1),PVT(:,2),60,PVT(:,3),'x')
            scatter(flipud(Pred+Presid),flipud(PVT(:,2)),80,flipud(PVT(:,3)),...
                'o','LineWidth',2)
            hold off;
        case 'normal'
            hold on;
            for(i=1:length(Tmod))
                iT = Tmod(i);
                iPmod = evalPressEos([],eosMod,Vmod,iT);
                h(i)=plot(iPmod,Vmod,'k-','LineWidth',2);
            end
            hold off;

            cmap=colormap();
            csamp = sampleColorMap(cmap,Tmod,opt.Tbnds);

            for(i=1:length(Tmod))
                set(h(i),'Color',csamp(i,:));
            end
            hold on;
            scatter(PVT(:,1),PVT(:,2),80,PVT(:,3),'o','LineWidth',2)
            hold off;
    end
    
    if(ismember(viewmode,{'normal','reduced'}))
        axes(opt.hax);
        caxis(opt.Tbnds);
        hcol=colorbar;
        
        set(hcol,'ytick',Tmod)
        ylim(opt.Vbnds)

        switch viewmode
            case 'reduced'
                xlabel('Reduced Pressure Isotherm [GPa]','FontSize',fontSz);
            case 'normal'
                xlabel('Pressure [GPa]','FontSize',fontSz);
        end
        ylabel(['Volume [' angChar '^3]'],'FontSize',fontSz);

        set(get(hcol,'ylabel'),'String','Temperature [K]','FontSize',fontSz)
        opt.hcolbar = hcol;
    end
    set(gca,'Box','on','FontSize',fontSz);
end
function optDefault = getDefaultOpt()
    optDefault.hax = [];
    optDefault.hcolbar = [];
    optDefault.Tbnds = [];
    optDefault.Vbnds = [];
end
