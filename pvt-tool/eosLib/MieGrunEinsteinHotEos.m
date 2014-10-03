% MieGrunEinsteinHotEos - calc. Einstein thermal contribution to P, K, Cv, etc
function [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunEinsteinHotEos(V,T,V0,T0,Natom,pHotEos,debyeDerivsFun)
    assert(false,'Mie-Gruneisen-Einstein Model not yet implemented');

    CvFac = pHotEos(end);

    % T in K
    % V in Ang^3
    % P in GPa
    R = 8.3144621/.6022*1e-3; 
    Cvm = 3*R*Natom*CvFac;

    % Gruneisen parameter is independent of Temp!
    [Tdeb,gam,dgamdV] = debyeDerivsFun(V,V0,pHotEos(1:end-1));
    gammaHot = gam;
    TdebyeHot= Tdeb;

    
   % EHot = Cvm.*(T - T0);
   % CvHot = Cvm*(4*debFunVal- 3*(Tdeb./T)./(exp(Tdeb./T)-1));

   % PHot = gam./V.*EHot;

   % dEHotdV = 3*PHot-3*Cvm*gam.*Tdeb./V.*...
   %     (1./(exp(Tdeb./T)-1) - 1./(exp(Tdeb./T0)-1));

   % KTHot = (gam./V-dgamdV).*EHot - gam.*dEHotdV;
end
