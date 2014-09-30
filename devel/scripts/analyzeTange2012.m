
angChar = char(197);

% Set reference temperature
T0 = 300;

% Analyze Tange (2012) pv data (using MgO press scale of Tange 2009)
datreadT12 = importdata('Tange2012pvdata.txt',' ',4);
%data_temp = cell2mat(textscan(fopen('test_file.txt'),'%f %f %f %f %f'));

run_colT12  = 1;
P_colT12    = 14;
Perr_colT12 = 15;
T_colT12    = 2;
Terr_colT12 = 3;
V_colT12    = 10;
Verr_colT12 = 11;
Vmark_colT12    = 12;
Vmarkerr_colT12 = 13;

pvdatT12 = initEosDataset('Tange2012',runID,...
    datreadT12.data(:,P_colT12), datreadT12.data(:,Perr_colT12),...
    datreadT12.data(:,V_colT12), datreadT12.data(:,Verr_colT12),...
    datreadT12.data(:,Vmark_colT12), datreadT12.data(:,Vmarkerr_colT12),...
    datreadT12.data(:,T_colT12), datreadT12.data(:,Terr_colT12),...
    0);

% Set 1 K errors for heated Sintered Diamond Multianvil
pvdatT12.Terr(pvdatT12.runID==1 & abs(pvdatT12.T/T0-1) > 3e-2) = 1;


pvdatT12.datsetName = 'Tange2012';
%Filter runID to 2 datatypes (multianvil and LHDAC)
pvdatT12.runID = round(datreadT12.data(:,run_colT12)/100);
pvdatT12.P        = datreadT12.data(:,P_colT12);
pvdatT12.Perr     = datreadT12.data(:,Perr_colT12);
pvdatT12.T        = datreadT12.data(:,T_colT12);
pvdatT12.Terr     = datreadT12.data(:,Terr_colT12);
pvdatT12.V        = datreadT12.data(:,V_colT12);
pvdatT12.Verr     = datreadT12.data(:,Verr_colT12);
pvdatT12.Vmark    = datreadT12.data(:,Vmark_colT12);
pvdatT12.Vmarkerr = datreadT12.data(:,Vmarkerr_colT12);
pvdatT12.PerrTot  = zeros(size(pvdatT12.P));
uniqIDT12 = unique(pvdatT12.runID);
pvdatT12.logErrModFac = zeros(length(uniqIDT12),3);

% Set 1 K errors for heated Sintered Diamond Multianvil
pvdatT12.Terr(pvdatT12.runID==1 & abs(pvdatT12.T/T0-1) > 3e-2) = 1;

datreadK09 = importdata('Katsura2009pvdata_corr.txt',' ',4);

run_colK09  = 1;
P_colK09    = 5;
Perr_colK09 = 6;
T_colK09    = 2;
Terr_colK09 = NaN;
V_colK09    = 7;
Verr_colK09 = 8;
Vmark_colK09    = 3;
Vmarkerr_colK09 = 4;

% NOTE: Need to scale Vmark and VmarkErr!
% This value was assumed for the volume of MgO in Katsura 2009
%   NOTE: there is some uncertainty of whether this was the precise value that
%   was used...
V0MgOK09 = 4.2112^3;
V0pvK09 = 4.7769*4.9298*6.8956;
pvdatK09.datsetName = 'Katsura2009';
pvdatK09.runID    = datreadK09.data(:,run_colK09);
pvdatK09.P        = datreadK09.data(:,P_colK09);
pvdatK09.Perr     = datreadK09.data(:,Perr_colK09);
pvdatK09.T        = datreadK09.data(:,T_colK09);
pvdatK09.Terr     = zeros(size(pvdatK09.T));
pvdatK09.V        = V0pvK09*datreadK09.data(:,V_colK09);
pvdatK09.Verr     = V0pvK09*datreadK09.data(:,Verr_colK09);
pvdatK09.Vmark    = V0MgOK09*datreadK09.data(:,Vmark_colK09);
pvdatK09.Vmarkerr = V0MgOK09*datreadK09.data(:,Vmarkerr_colK09);
pvdatK09.PerrTot  = zeros(size(pvdatK09.P));
uniqIDK09 = unique(pvdatK09.runID);
pvdatK09.logErrModFac = zeros(length(uniqIDK09),3);

% Combine all data into list
pvdatAll = [pvdatT12 pvdatK09];



% Construct EOS objects for initial guess and priors
pveos_init.T0  = T0;
pveos_init.coldFun = @eosVIN;
pveos_init.MGDFun  = @calcPthmMGDpowlaw;
pveos_init.markcoldFun = @eosVIN;
pveos_init.markMGDFun  = @calcPthmMGDpowlaw;
%                      V0     K0    KP0  TD gam0  q  Nat
pveos_init.eos = [162.373 258.4 4.10 940 1.55 1.1 4*5];
pveos_init.eoscov = diag([0 1.7 0.07 140 0.09 0.3 0].^2);
% Reported values from Tange 2009 (Fit2 - Vinet)
pveos_init.eosmark = [74.698 160.61 4.454 762 1.451 0.86 4*2];

% Set priors based on previous measurements and physical limits
pveos_prior = pveos_init;
pveos_prior.eos = [162.5 258.4 4.10 1100 1.0 1.0 4*5];
pveos_prior.eoscov = diag([0.2 Inf Inf 0 1.0 1.0 0].^2);

% Alternate prior that allows fitting of Debye Temp
pveos_priorDeb = pveos_prior;
pveos_priorDeb.eos(4) = 1000;
pveos_priorDeb.eoscov(4,4) = 400.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fit EOS for Tange 2012 + Katsura 2009 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvdatAllR = pvdatAll;

doRobustFit=false;
updatePscl = true;
pvdatAll  = updateTotErrstruc(pveos_init,pvdatAll,doRobustFit,updatePscl);
[pveosAll_fit jointDatStruc_fit] = fitEOSstruc(pveos_init,pveos_prior,pvdatAll,doRobustFit);

doRobustFit=true;
updatePscl = true;
pvdatAllR  = updateTotErrstruc(pveos_init,pvdatAllR,doRobustFit,updatePscl);
[pveosAll_fitR jointDatStruc_fitR] = fitEOSstruc(pveos_init,pveos_prior,pvdatAllR,doRobustFit);
%pvdatAllR  = updateTotErrstruc(pveosT12_init,pvdatAll,doRobustFit,updatePscl);


% Recalculate with Debye Temp as free parameter
doRobustFit=true;
updatePscl = true;
pvdatAllR  = updateTotErrstruc(pveos_init,pvdatAllR,doRobustFit,updatePscl);
[pveosAll_fitDebR jointDatStruc_fitDebR] = fitEOSstruc(pveos_init,pveos_priorDeb,pvdatAllR,doRobustFit);

pveosAll_fitDebR.eos([1:6])
sqrt(diag(pveosAll_fitDebR.eoscov(1:6,1:6)))'
pvdatAllR.PerrTot

V0=pveosT12_init.eos(1);
V300 = V0*linspace(1,.75,300);;
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pveosT12_init.eos(1:3));
[Pmod300F,dEmod300F,Kmod300F,KPmod300F] = eosVIN(V300,pveosAll_fit.eos(1:3));
clf;hold on;
scatter(pvdatAll(1).P,pvdatAll(1).V/V0,50,pvdatAll(1).T,'o')
scatter(pvdatAll(2).P,pvdatAll(2).V/V0,50,pvdatAll(2).T,'x')
plot(0,1,'bx')
plot(Pmod300,V300/V0,'b:');
plot(Pmod300F,V300/V0,'b-');
hold off;
colorbar;
caxis([300 3000])
[pveosAll_fit jointDatStruc_fit] = fitEOSstruc(pveosT12_init,pveosT12_prior,pvdatAll,doRobustFit);


doRobustFit=true;
updatePscl = true;
pvdatAllR  = updateTotErrstruc(pveosT12_init,pvdatAll,doRobustFit,updatePscl);
[pveosAll_fitR jointDatStruc_fitR] = fitEOSstruc(pveosT12_init,pveosT12_prior,pvdatAllR,doRobustFit);

pvdatAll  = updateTotErrstruc(pveosT12_init,pvdatAll,doRobustFit,updatePscl);

pvdatT12  = updateTotErrstruc(pveosT12_init,pvdatT12,doRobustFit,updatePscl);

[pveosT12_fit residT12] = fitEOSstruc(pveosT12_init,pveosT12_prior,pvdatT12,doRobustFit);
pvdatT12  = updateTotErrstruc(pveosT12_fit,pvdatT12,doRobustFit,updatePscl);
[pveosT12_fit residT12] = fitEOSstruc(pveosT12_fit,pveosT12_prior,pvdatT12,doRobustFit);

doRobustFit=true;
pvdatT12R = pvdatT12;
pveosT12R_init = pveosT12_init;
pveosT12R_prior = pveosT12_prior;
pvdatT12R  = updateTotErrstruc(pveosT12R_init,pvdatT12R,doRobustFit,updatePscl);
[pveosT12R_fit residT12R] = fitEOSstruc(pveosT12R_init,pveosT12R_prior,pvdatT12R,doRobustFit);
pvdatT12R  = updateTotErrstruc(pveosT12R_fit,pvdatT12R,doRobustFit,updatePscl);
[pveosT12R_fit residT12R] = fitEOSstruc(pveosT12R_fit,pveosT12R_prior,pvdatT12R,doRobustFit);

% Compare robust and normal fits

% Cold params
pveosT12_fit.eos(1:3)
pveosT12R_fit.eos(1:3)
sqrt(diag(pveosT12_fit.eoscov(1:3,1:3)))'
sqrt(diag(pveosT12R_fit.eoscov(1:3,1:3)))'

% Hot params
pveosT12_fit.eos(5:6)
pveosT12R_fit.eos(5:6)
sqrt(diag(pveosT12_fit.eoscov(5:6,5:6)))'
sqrt(diag(pveosT12R_fit.eoscov(5:6,5:6)))'

% They agree to within uncertainties


% Recalculate with Debye Temp as free parameter
pveosT12R_priorDeb = pveosT12R_prior;
pveosT12R_priorDeb.eos(4) = 1000;
pveosT12R_priorDeb.eoscov(4,4) = 400.^2;

doRobustFit = true;
[pveosT12R_fitDeb residT12RDeb] = fitEOSstruc(pveosT12R_fit,pveosT12R_priorDeb,pvdatT12R,doRobustFit);
% Need to update var again?

% Compare all fits

% Cold params
pveosT12_fit.eos(1:3)
pveosT12R_fit.eos(1:3)
pveosT12R_fitDeb.eos(1:3)
sqrt(diag(pveosT12_fit.eoscov(1:3,1:3)))'
sqrt(diag(pveosT12R_fit.eoscov(1:3,1:3)))'
sqrt(diag(pveosT12R_fitDeb.eoscov(1:3,1:3)))'

% Hot params
pveosT12_fit.eos(5:6)
pveosT12R_fit.eos(5:6)
pveosT12R_fitDeb.eos(5:6)
sqrt(diag(pveosT12_fit.eoscov(5:6,5:6)))'
sqrt(diag(pveosT12R_fit.eoscov(5:6,5:6)))'
sqrt(diag(pveosT12R_fitDeb.eoscov(5:6,5:6)))'

pveosT12R_fitDeb.eos(4)
sqrt(pveosT12R_fitDeb.eoscov(4,4))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fit EOS for Katsura 2009 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update Tot Errors based on initial guess
%
doRobustFit = false;
updatePscl = true;
pvdatK09 = updateTotErrstruc(pveosK09_init,pvdatK09,doRobustFit,updatePscl);

doRobustFit=false;
[pveosK09_fit residK09] = fitEOSstruc(pveosK09_init,pveosK09_prior,pvdatK09,doRobustFit);

%Re-update Tot Errors using fit
pvdatT12 = updateTotErrstruc(pveosT12_init,pvdatT12);

doRobustFit=false;
[pveosT12_fit residT12] = fitEOSstruc(pveosT12_init,pveosT12_prior,pvdatT12,doRobustFit);

doRobustFit=true;
[pveosT12_fitR residT12R] = fitEOSstruc(pveosT12_init,pveosT12_prior,pvdatT12,doRobustFit);

% Recalculate with Debye Temp as free parameter
pveosT12_priorDeb = pveosT12_prior;
pveosT12_priorDeb.eos(4) = 1000;
pveosT12_priorDeb.eoscov(4,4) = 200.^2;

doRobustFit=false;
[pveosT12_fitDeb residT12Deb] = fitEOSstruc(pveosT12_fit,pveosT12_priorDeb,pvdatT12,doRobustFit);

% Compare fits
sqrt(diag(pveosT12_fitDeb.eoscov))
sqrt(diag(pveosT12_priorDeb.eoscov))
sqrt(diag(pveosT12_fit.eoscov))
pveosT12_fitDeb
pveosT12_fit
pveosT12_fitR




% First update Tot Errors based on initial guess
pvdatK09 = updateTotErrstruc(pveosK09_init,pvdatK09);

doRobustFit=false;
[pveosK09_fit residK09] = fitEOSstruc(pveosK09_init,pveosK09_prior,pvdatK09,doRobustFit);

pveosT12_fit     = fitEOSstruc(pveosT12_fitCold,pveosT12_fitCold,pvdatT12,...
    find(pvdatT12.T>300),doRobustFit);


pveosT12_pri.eos = [
pveosT12_pri.eoscov
pveosT12_fitInit = fitColdEOSstruc(pveosT12_pri,pveosT12_pri,pvdatT12);
pveosT12_fit = estTotErrColdstruc(pveosT12_fitInit,pveosT12_pri,pvdatT12);


pveosK09.eos
pveosK09.eoscov
%pveosK09_fit = fitColdEOSstruc(pveosT12_pri,pveosT12_pri,pvdatT12);
pveosK09 = pveosK09_pri;




% NOTE: This is the "Fit2" option from Tange 2009
%   -> This has a standard MGD EOS form (a fixed to 1)
peosMgOT09 = [74.698 160.61 4.454 762 1.451 0.86 2*4];
eosHotFun = @(V,T,peos)(eosVIN(V,peos(1:3))+...
    calcPthmMGDpowlaw(V,T,peos(1),300,peos(4:end)));
[Pmark,dPdVmark,dPdTmark] = evalEosPVT(Vmark,T,peosMgOT09,eosHotFun);


runtyp=floor(runID/100);

scatter(P,V,50,T,'o','filled')
colorbar;

ind300=find(T==300);
indHot=find(T>300);

runtyp300 = runtyp(ind300);
runID300 = runID(ind300);
P300 = P(ind300);
V300 = V(ind300);
Vmark300 = Vmark(ind300);
Perr300 = Perr(ind300);
Verr300 = Verr(ind300);
Vmarkerr300 = Vmarkerr(ind300);

%PerrTot300 = sqrt(Perr300.^2 +  (dPdVmark(ind300).*Vmarkerr300).^2);
PerrTot300 = sqrt((dPdVmark(ind300).*Vmarkerr300).^2);

doRobustFit=false;
prior = [162.373 258 4.1];
prior = [162.5 260 4.1];
priorcov = Inf*ones(size(prior));
priorcov(1) = (1e-6*prior(1))^2;
pinit = prior;
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);
plot(P300,V300,'bo',Pmod300,V300,'rx')

%%%%%%%%%%%%%%%%%%
%  Use Wolf priors, calc propagated error
%     -> Tange reported errors that are decent to order of mang (300K)
%        but they still underrepresent errors slightly
%%%%%%%%%%%%%%%%%%
%prior = [162.5 260 4.1];
priorcov = Inf*ones(size(prior));
%priorcov(1) = [.2^2];
priorcov(1) = [1e-4^2];
doRobustFit=false;
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);

plot(P300,V300,'bo',Pmod300,V300,'rx')
plot(P300,Kmod300,'bo')

% Add error contribution from error in volume
%PerrTot300 = sqrt(Perr300.^2 + abs(Kmod300.*Verr300./V300).^2 + ...
%    (dPdVmark(ind300).*Vmarkerr300).^2);
PerrTot300 = sqrt(abs(Kmod300.*Verr300./V300).^2 + ...
    (dPdVmark(ind300).*Vmarkerr300).^2);
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);
1.4826*median(abs((Pmod300-P300)./PerrTot300))

% Separate out 2 diff datatypes
1.4826*median(abs((Pmod300(runtyp300==1)-P300(runtyp300==1))./PerrTot300(runtyp300==1)))
1.4826*median(abs((Pmod300(runtyp300==2)-P300(runtyp300==2))./PerrTot300(runtyp300==2)))

%sqrt(mean(((Pmod300(runtyp300==1)-P300(runtyp300==1))./PerrTot300(runtyp300==1)).^2))
%sqrt(mean(((Pmod300(runtyp300==2)-P300(runtyp300==2))./PerrTot300(runtyp300==2)).^2))

priorcov(1) = [.2^2];
prior = pcoldFit;
pinit = pcoldFit;
prior(1) = 162.5;
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);

1.4826*median(abs((Pmod300-P300)./PerrTot300))

% Separate out 2 diff datatypes
1.4826*median(abs((Pmod300(runtyp300==1)-P300(runtyp300==1))./PerrTot300(runtyp300==1)))
1.4826*median(abs((Pmod300(runtyp300==2)-P300(runtyp300==2))./PerrTot300(runtyp300==2)))


Ndraw = 100;
pdraw = mvnrnd(pcoldFit,pcoldCov,Ndraw);
clf;
ploterr(P300,V300,Perr300,Verr300,'bo','hhxy',0)
Varr = linspace(min(V300),max(V300),60);
hold on;
for(i=1:Ndraw);
    plot(eosVIN(Varr,pdraw(i,:)),Varr,'k-');
end
ploterr(P300,V300,Perr300,Verr300,'bo','hhxy',0)
ploterr(P300,V300,PerrTot300,[],'ro','hhxy',0)
hold off;

%%%%%%%%%%%%%%%%%%%%%%%
%  Inflate errors
%%%%%%%%%%%%%%%%%%%%%%%
logErrVfacArr = linspace(-6,-1,100);
PscatMad1= zeros(size(logErrVfacArr));
PscatStd1= zeros(size(logErrVfacArr));
PscatMad2= zeros(size(logErrVfacArr));
PscatStd2= zeros(size(logErrVfacArr));
for(i=1:length(logErrVfacArr))
    iVfac = exp10(logErrVfacArr(i));
    %iPerrTot300 = sqrt(Perr300.^2 +(dPdVmark(ind300).*(Vmarkerr300+iVfac*Vmark300)).^2 + ...
    %    (Kmod300.*(Verr300./V300)).^2 );
    iPerrTot300 = sqrt((dPdVmark(ind300).*(Vmarkerr300+iVfac*Vmark300)).^2 + ...
        (Kmod300.*(Verr300./V300)).^2 );
    PscatMad1(i)=1.4826*median(abs((Pmod300(runtyp300==1)-P300(runtyp300==1))./iPerrTot300(runtyp300==1)));
    PscatStd1(i)=sqrt(mean(((Pmod300(runtyp300==1)-P300(runtyp300==1))./iPerrTot300(runtyp300==1)).^2));
    PscatMad2(i)=1.4826*median(abs((Pmod300(runtyp300==2)-P300(runtyp300==2))./iPerrTot300(runtyp300==2)));
    PscatStd2(i)=sqrt(mean(((Pmod300(runtyp300==2)-P300(runtyp300==2))./iPerrTot300(runtyp300==2)).^2));
end
%plot(logErrVfacArr,log10(PscatMad),'r-',logErrVfacArr,log10(PscatStd),'k-')
plot(logErrVfacArr,log10(PscatMad1),'k-',logErrVfacArr,log10(PscatStd1),'k--',...
    logErrVfacArr,log10(PscatMad2),'r-',logErrVfacArr,log10(PscatStd2),'r--');

%errVfac =exp10(interp1(log10(PscatStd),logErrVfacArr,0));
%errVfacR=exp10(interp1(log10(PscatMad),logErrVfacArr,0));
errVfacR1=exp10(interp1(log10(PscatMad1),logErrVfacArr,0));
errVfacR2=exp10(interp1(log10(PscatMad2),logErrVfacArr,0));

%PerrTot300 = sqrt(Perr300.^2 +(dPdVmark(ind300).*(Vmarkerr300+errVfacR*Vmark300)).^2 + ...
%    (Kmod300.*(Verr300./V300)).^2 );
PerrTot300(runtyp300==1) = sqrt((dPdVmark(ind300(runtyp300==1)).*(Vmarkerr300(runtyp300==1)+errVfacR1*Vmark300(runtyp300==1))).^2 + ...
    (Kmod300(runtyp300==1).*(Verr300(runtyp300==1)./V300(runtyp300==1))).^2 );
PerrTot300(runtyp300==2) = sqrt((dPdVmark(ind300(runtyp300==2)).*(Vmarkerr300(runtyp300==2)+errVfacR2*Vmark300(runtyp300==2))).^2 + ...
    (Kmod300(runtyp300==2).*(Verr300(runtyp300==2)./V300(runtyp300==2))).^2 );
%PerrTot300 = sqrt(Perr300.^2 + abs(Kmod300.*(Verr300./V300+errVfacR)).^2);
%1.4826*median(abs((Pmod300-P300)./PerrTot300))
%sqrt(mean(((Pmod300-P300)./PerrTot300).^2))

% Separate out 2 diff datatypes
1.4826*median(abs((Pmod300(runtyp300==1)-P300(runtyp300==1))./PerrTot300(runtyp300==1)))
1.4826*median(abs((Pmod300(runtyp300==2)-P300(runtyp300==2))./PerrTot300(runtyp300==2)))

%prior = [162.5 260 4.1];
pinit = prior;
%prior = [162.373 260 4.1];
priorcov = Inf*ones(size(prior));
priorcov(1) = [.2^2];
%priorcov(1) = [1e-4^2];
doRobustFit=false;
%doRobustFit=true;
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);

Ndraw = 300;
Nsamp = 100;
pdraw = mvnrnd(pcoldFit,pcoldCov,Ndraw);
Vmod = linspace(min(V300),max(V300),Nsamp);
Pmoddraw = zeros(Ndraw,length(Vmod));
for(i=1:Ndraw);
    Pmoddraw(i,:)=eosVIN(Vmod,pdraw(i,:));
end
PmodBnds = quantile(Pmoddraw,normcdf([-1 1]));
clf;
hBnds=fill([PmodBnds(1,:) fliplr(PmodBnds(2,:))],[Vmod fliplr(Vmod)],[.8 .8 .95]);
set(hBnds,'EdgeColor','none');
hold on;
herr=ploterr(P300,V300,PerrTot300,[],'bo','hhxy',0)
set(herr,'MarkerEdgeColor','none','LineWidth',2)
hold off;
xlabel('Pressure [GPa]')
ylabel(['Volume [' angChar '^3]'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Fit Hot EOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%indHot=find(T>300);
PHot = P(indHot);
VHot = V(indHot);
THot = T(indHot);
PerrHot = Perr(indHot);
VerrHot = Verr(indHot);
TerrHot = Terr(indHot);
VmarkHot    = Vmark(indHot);
VmarkerrHot = Vmarkerr(indHot);

runtypHot = runtyp(indHot);
runIDHot = runID(indHot);

%[Pmark,dPdVmark,dPdTmark] = evalEosPVT(Vmark,T,peosMgOT09,eosHotFun);
%  DO NOT INCLUDE PERRHOT... it is double counting errors
%PerrTotHot = sqrt(PerrHot.^2 +  (dPdVmark(indHot).*VmarkerrHot).^2 + ...
%    (dPdTmark(indHot).*TerrHot).^2);
PerrTotHot = sqrt((dPdVmark(indHot).*VmarkerrHot).^2 + ...
    (dPdTmark(indHot).*TerrHot).^2);

priorHot    = [1100    1.0 1.0 20];
%priorHot    = [ 940    1.0 1.0 20];
priorWidHot = [0       1 1 0];

priorAll = [pcoldFit priorHot];
priorcovAllinit      = zeros(3+length(priorHot));
priorcovAllinit(4,4) = priorWidHot(1)^2;
priorcovAllinit(5,5) = priorWidHot(2)^2;
priorcovAllinit(6,6) = priorWidHot(3)^2;
priorcovAllinit(7,7) = priorWidHot(4)^2;


doRobustFit = false;
[pallFitInit,pallCov,residHotInit] = fitHotEosPVT(priorAll,priorAll,...
    priorcovAllinit,PHot,VHot,THot,PerrTotHot,doRobustFit,@eosVIN);

%eosHotFun = @(V,T,peos)(eosVIN(V,peos(1:3))+...
%    calcPthmMGDpowlaw(V,T,peos(1),300,peos(4:end)));

[PHot-eosHotFun(VHot,THot,pallFitInit),residHotInit]

[PmodHotInit,dPdVmodHotInit,dPdTmodHotInit] = evalEosPVT(VHot,THot,...
    pallFitInit,eosHotFun);

% Evaluate propagated errors
%PerrTotHot = sqrt(PerrHot.^2 +  (dPdVmark(indHot).*VmarkerrHot).^2 + ...
%    ((dPdTmark(indHot)-dPdTmodHotInit).*TerrHot).^2 + ...
%    (dPdVmodHotInit.*VerrHot).^2);
PerrTotHot = sqrt((dPdVmark(indHot).*VmarkerrHot).^2 + ...
    ((dPdTmark(indHot)-dPdTmodHotInit).*TerrHot).^2 + ...
    (dPdVmodHotInit.*VerrHot).^2);

doRobustFit = false;
doRobustFit = true;
priorcovAll = priorcovAllinit;
priorcovAll(1:3,1:3)     = pcoldCov;
[pallFit,pallCov,residHot] = fitHotEosPVT(priorAll,priorAll,...
    priorcovAll,PHot,VHot,THot,PerrTotHot,doRobustFit,@eosVIN);

[PmodHot,dPdVmodHot,dPdTmodHot] = evalEosPVT(VHot,THot,pallFit,eosHotFun);
PerrTotHot = sqrt((dPdVmark(indHot).*VmarkerrHot).^2 + ...
    ((dPdTmark(indHot)-dPdTmodHot).*TerrHot).^2 + ...
    (dPdVmodHot.*VerrHot).^2);

[PHot-PmodHot,residHot]
1.4826*median(abs((PmodHot-PHot)./PerrTotHot))

1.4826*median(abs((PmodHot(runtypHot==1)-PHot(runtypHot==1))./PerrTotHot(runtypHot==1)))
1.4826*median(abs((PmodHot(runtypHot==2)-PHot(runtypHot==2))./PerrTotHot(runtypHot==2)))

% Propagated errors are already a bit high

sqrt(mean(((PmodHot-PHot)./PerrTotHot).^2))

residHot./PerrTotHot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   No need to inflate errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ndraw = 50;
Nsamp = 45;
palldraw = mvnrnd(pallFit,pallCov,Ndraw);
Vmod = linspace(min(V),max(V),Nsamp);
Tmod = 2500;
Pmoddraw = zeros(Ndraw,length(Vmod));
for(i=1:Ndraw);
    [iPmoddraw,dPdVmoddraw,dPdTmoddraw] = evalEosPVT(Vmod,Tmod,palldraw(i,:),eosHotFun);
    Pmoddraw(i,:) = iPmoddraw;
    %Pmoddraw(i,:)=eosVIN(Vmod,pdraw(i,:));
end
PmodBnds = quantile(Pmoddraw,normcdf([-1 1]));

clf;
%hBnds=fill([PmodBnds(1,:) fliplr(PmodBnds(2,:))],[Vmod fliplr(Vmod)],[.8 .8 .95]);
hBnds=fill([PmodBnds(1,:) fliplr(PmodBnds(2,:))],[Vmod fliplr(Vmod)],[.8 .95 .8]);
hBnds=fill([PmodBnds(1,:) fliplr(PmodBnds(2,:))],[Vmod fliplr(Vmod)],[.95 .8 .8]);
set(hBnds,'EdgeColor','none');
hold on;
herr=ploterr(P300,V300,PerrTot300,[],'b.','hhxy',0)
herrHot=ploterr(PHot,VHot,PerrTotHot,[],'r.','hhxy',0)
set(herr,'MarkerEdgeColor','none','LineWidth',2)
set(herrHot,'MarkerEdgeColor','none','LineWidth',2)
hold off;
xlabel('Pressure [GPa]')
ylabel(['Volume [' angChar '^3]'])

hold on;
hscat=scatter(P,V,80,T,'o')
set(hscat,'LineWidth',1.5)
colorbar;
hold off;

pallFit(1:3)
pallFit(5:6)



scatter(PmodHot,VHot,80,THot,'o')
hold on;
scatter(PHot,VHot,80,THot,'x')
hold off;



scatter(PmodHotInit,VHot,80,THot,'o')
hold on;
scatter(PHot,VHot,80,THot,'x')
hold off;

[PvarMeas,PvarVmark] = initEmpiricalTotErr(Vmark,dPdVmark,...
    VsampErr,dPdVsamp, Terr,dPdTsamp,dPdTmark)

[PvarMeasHotW14,PvarVmarkHotW14] = initEmpiricalTotErr(VmarkMeasW14(indHotW14,1),...
    PderivW14(indHotW14,1),VerrW14(indHotW14,1),dPdVmodHotInit,...
    VTerrW14(indHotW14,2),dPdTmodHot87init,PderivW14(indHotW14,2));



[PmodHot87init,dPdVmodHot87init,dPdTmodHot87init] = evalEosPVT(PVTW14(indHotW14,2),...
    PVTW14(indHotW14,3),pall87init,eosHotFun);

[PvarMeasHot,PvarVmarkHot] = initEmpiricalTotErr(VmarkMeas(indHot,1),...
    Pderiv(indHot,1),VTerr(indHot,1),dPdVmodHotinit,...
    VTerr(indHot,2),dPdTmodHotinit,Pderiv(indHot,2));


[PerrTotHot,logVerrFacHot] = getEmpiricalTotErr(residHotInit,...
    PvarMeasHot,PvarVmarkHot);

pallFit(1:3)
pallFit(5:6)

%priorcovAllinit(1:3,1:3) = pcoldCov;
pallErr = sqrt(diag(pallCov));
pallCorr = pallCov./(pallErr(:)*pallErr(:)');


scatter(P,V,50,T,'o','filled')
colorbar;




[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);
plot(P300,V300,'bo',Pmod300,V300,'rx')
plot(P300,Kmod300,'bo')
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);
1.4826*median(abs((Pmod300-P300)./PerrTot300))
priorcov(1) = [.2^2];
prior = pcoldFit;
pinit = pcoldFit;
prior(1) = 162.5;
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);
[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);
1.4826*median(abs((Pmod300-P300)./PerrTot300))
sqrt(mean(((Pmod300-P300)./PerrTot300).^2))




%doRobustFit=true;
%[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P300,V300,PerrTot300,doRobustFit,@eosVIN);
%[Pmod300,dEmod300,Kmod300,KPmod300] = eosVIN(V300,pcoldFit);
%1.4826*median(abs((Pmod300-P300)./PerrTot300))
%sqrt(mean(((Pmod300-P300)./PerrTot300).^2))

sqrt(diag(pcoldCov))
sqrt(diag(pcoldCov))


[pall87init,pallcov87init,residHot87init] = fitHotEosPVT(priorAll87,priorAll87,...
    priorcovAll87init,PVTW14(indHotW14,1),PVTW14(indHotW14,2),PVTW14(indHotW14,3),...
    [],false,@eosVIN);

plot(P300,V300,'bo',Pmod300,V300,'rx')


plot(P300,V300,'bo',Pmod300,V300,'rx')
plot(P300,Kmod300,'bo')

[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P,V,Perr,doRobustFit,@eosVIN);
[Pmod,dEmod,Kmod,KPmod] = eosVIN(V,pcoldFit);

plot(P,V
resid./Perr

%[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pcoldFit,prior,priorcov,P,V,Perr,doRobustFit,@eosVIN);
[pcoldFit,pcoldCov,resid,nLogLkMin] = fitColdEosPV(pinit,prior,priorcov,P,V,Perr,doRobustFit,@eosVIN);

ploterr(P(ind300),V(ind300),Perr(ind300),Verr(ind300),'ko')
fit

