dataDir = '/Users/aswolf/Documents/code/MATLAB/pvt-tool/devel/data';
figDir = '/Users/aswolf/Documents/code/MATLAB/pvt-tool/devel/figs';


% Fit MgPv (Tange2012) and MgFePv (Wolf2014) datasets
lsqMgPv  =runPVTtool(fullfile(dataDir,'fit_Tange2012_lsq_fixdeb.in'));
lsqMgPvNothmgrp =runPVTtool(fullfile(dataDir,'fit_Tange2012_lsq_fixdeb_nothmgrp.in'));
lsqMgFePv=runPVTtool(fullfile(dataDir,'fit_Wolf2014_lsq_fixdeb.in'));
lsqMgFePv1=runPVTtool(fullfile(dataDir,'fit_test_grp1.in'));

fig=1;
figure(fig);
clf;
viewPVTFit(lsqMgFePv,'normal')
exportfigmixedrender(fig,fullfile(figDir,'Wolf14-PVT-fit.eps'))
viewPVTFit(lsqMgFePv,'reduced')
exportfigmixedrender(fig,fullfile(figDir,'Wolf14-reduced-fit.eps'))
viewPVTFit(lsqMgFePv,'hist')
exportfigmixedrender(fig,fullfile(figDir,'Wolf14-hist-resid.eps'))

viewPVTFit(lsqMgPv,'normal')
exportfigmixedrender(fig,fullfile(figDir,'Tange12-PVT-fit.eps'))
viewPVTFit(lsqMgPv,'reduced')
exportfigmixedrender(fig,fullfile(figDir,'Tange12-reduced-fit.eps'))
viewPVTFit(lsqMgPv,'hist')
exportfigmixedrender(fig,fullfile(figDir,'Tange12-hist-resid.eps'))




%LOAD PREM data
PREMload=importdata('~/Documents/data/PREM1s.txt');
PREM = PREMload.data;
lowMantleInd = find(PREM(:,1)==3|PREM(:,1)==4);
[tr,indUniq] = unique(PREM(lowMantleInd,4));
Pprem=PREM(lowMantleInd(indUniq),4);
rhoprem=PREM(lowMantleInd(indUniq),7);
zprem=PREM(lowMantleInd(indUniq),3);
Kprem = PREM(lowMantleInd(indUniq),8);
vcprem=sqrt(PREM(lowMantleInd(indUniq),8)./PREM(lowMantleInd(indUniq),7));
Pcmb = 135.8;
vcpremcmb=interp1(Pprem,vcprem,Pcmb,'spline');
rhopremcmb=interp1(Pprem,rhoprem,Pcmb,'spline');
Kpremcmb = interp1(Pprem,Kprem,Pcmb,'spline');

PpremHR   = linspace(Pprem(1),Pcmb,100);
rhopremHR = interp1(Pprem,rhoprem,PpremHR,'spline','extrap');
plot(Pprem,rhoprem,'k-o',Pprem,rhopremcmb*(1+(Pprem-Pcmb)./Kpremcmb),'r--',PpremHR,rhopremHR,'gx')


% Examine CMB properties of each composition
Pbnds = [30, 130];
Tbnds = [300, 3000];
eosModEndmem = [lsqMgPv.sampEosFit, lsqMgFePv.sampEosFit];
Xendmem = [0, 0.13];

TCMB=2450;
PCMB=135.8;
V0 = 160;
evalPMg = @(V,T)(evalPressEos([],eosModEndmem(1),V,T));
VCMB_Mg = fzero(@(V)(evalPMg(V,TCMB)-PCMB),V0);
evalPMgFe = @(V,T)(evalPressEos([],eosModEndmem(2),V,T));
VCMB_MgFe = fzero(@(V)(evalPMgFe(V,TCMB)-PCMB),V0);
[PCMB_Mg,PderivsCMB_Mg,KTCMB_Mg,CvCMB_Mg,gamCMB_Mg,thmExpCMB_Mg] = evalPMg(VCMB_Mg,TCMB);
KSCMB_Mg = KTCMB_Mg*(1+thmExpCMB_Mg*gamCMB_Mg*TCMB);
[PCMB_MgFe,PderivsCMB_MgFe,KTCMB_MgFe,CvCMB_MgFe,gamCMB_MgFe,thmExpCMB_MgFe] = evalPMgFe(VCMB_MgFe,TCMB);
KSCMB_MgFe = KTCMB_MgFe*(1+thmExpCMB_MgFe*gamCMB_MgFe*TCMB);
100*(KSCMB_MgFe./KSCMB_Mg-1)
100*(KTCMB_MgFe./KTCMB_Mg-1)
100*(thmExpCMB_MgFe./thmExpCMB_Mg-1)
100*(KSCMB_MgFe./Kpremcmb-1)
100*(KSCMB_Mg./Kpremcmb-1)
100*(KTCMB_MgFe./Kpremcmb-1)
100*(KTCMB_Mg./Kpremcmb-1)
Tfoot = TCMB;
dlogV = .01;
[PadMg,VadMg,TadMg,KadMg,gamadMg,thmExpadMg]=calcAdiabat(Tfoot,PCMB,30,dlogV,eosModEndmem(1));
[PadMgFe,VadMgFe,TadMgFe,KadMgFe,gamadMgFe,thmExpadMgFe]=calcAdiabat(Tfoot,PCMB,30,dlogV,eosModEndmem(2));
plot(PadMg,VadMg,'k-',PadMgFe,VadMgFe,'r-')
plot(PadMg,KadMg,'k-',PadMgFe,KadMgFe,'r-')

Xideal = [0:.01:0.2];
eosModIdeal = calcIdealLatticeMix(Xideal, Xendmem, eosModEndmem,Pbnds,Tbnds);

% pEosIdealFit = reshape([eosModIdeal.pEos],[],length(eosModIdeal))';
% plot(Xideal,pEosIdealFit(:,1),'k-')
% plot(Xideal,pEosIdealFit(:,2),'k-')
% plot(Xideal,pEosIdealFit(:,3),'k-')
% plot(Xideal,pEosIdealFit(:,5),'k-')
% plot(Xideal,pEosIdealFit(:,6),'k-')

m100 = 24.31+28.09+3*16;
m0 = 55.85+28.09+3*16;
m87 = .87*24.31+.13*55.85+28.09+3*16;
mendmem = [m100, m87];
Z = 4;
wtendmem = (Xideal-Xendmem(1))./(Xendmem(2)-Xendmem(1));
mmix = Z*((1-wtendmem).*mendmem(1) + wtendmem.*mendmem(2))/.6022;

TCMB=2450;
PCMB=135.8;
Pfoot=PCMB;
Pstop = 30;
dlogV = .01;
TpileExv = [0:100:1200];
Pad = linspace(30,PCMB,100);
VadM      = zeros(length(Xideal), length(TpileExv), length(Pad));
rhoadM    = zeros(length(Xideal), length(TpileExv), length(Pad));
drhoM    = zeros(length(Xideal), length(TpileExv), length(Pad));
dKadM    = zeros(length(Xideal), length(TpileExv), length(Pad));
TadM      = zeros(length(Xideal), length(TpileExv), length(Pad));
KadM      = zeros(length(Xideal), length(TpileExv), length(Pad));
thmExpadM = zeros(length(Xideal), length(TpileExv), length(Pad));
% clf; hold on;
rhopremv = interp1(Pprem,rhoprem,Pad,'spline','extrap');
Kpremv = interp1(Pprem,Kprem,Pad,'spline','extrap');
for(i=1:length(Xideal))
    i/length(Xideal)
    iX = Xideal(i);
    ieosMod = eosModIdeal(i);
    for(j=1:length(TpileExv))
        jTpileEx = TpileExv(j);
        jTfoot=TCMB+jTpileEx;
        [ijPad,ijVad,ijTad,ijKad,ijgamad,ijthmExpad]=...
            calcAdiabat(jTfoot,Pfoot,Pstop,dlogV,ieosMod);
        ijrhoad = mmix(i)./ijVad;
        VadM(i,j,:)      = interp1(ijPad,ijVad,Pad,'spline','extrap');
        rhoadM(i,j,:)    = interp1(ijPad,ijrhoad,Pad,'spline','extrap');
        drhoM(i,j,:)     = squeeze(rhoadM(i,j,:))./rhopremv' - 1;
        dKadM(i,j,:)     =interp1(ijPad,ijKad,Pad,'spline','extrap')./Kpremv-1;
        TadM(i,j,:)      = interp1(ijPad,ijTad,Pad,'spline','extrap');
        KadM(i,j,:)      = interp1(ijPad,ijKad,Pad,'spline','extrap');
        thmExpadM(i,j,:) = interp1(ijPad,ijthmExpad,Pad,'spline','extrap');
        %plot(Pad,ijrhoad,'k-',Pprem,rhoprem,'r-')
        % plot(Pad,100*(ijrhoad./interp1(Pprem,rhoprem,Pad,'spline','extrap')-1),'k-');
        % pause(.01)
    end
end



drhoCMB = squeeze(drhoM(:,:,end));
dKCMB = squeeze(dKadM(:,:,end));
contour(Xideal, TpileExv,100*dKCMB',[-10:10]);colorbar;axis xy;
contour(Xideal, TpileExv,100*drhoCMB',[-10:10]);colorbar;axis xy;
contour(Xideal, TpileExv,1e5*thmExpadM(:,:,end)');colorbar;axis xy;

% 40 to 68

Pnb = zeros(length(Xideal),length(TpileExv));
for(i=1:length(Xideal))
    for(j=1:length(TpileExv))
        %plot(Pad,squeeze(rhoadM(i,j,:)),'ko-',Pad,rhopremv,'ro-')
        Pnb(i,j) = interp1(Pad,drhoM(i,j,:),
        plot(Pad,100*squeeze(drhoM(i,j,:)),'ko-')
        pause(.01)
    end
end



plot(Pad,mmix(ind)./Vad,'k-')

plot(Pad,Tad,'k-');
plot(Pad,1./Vad,'k-');
