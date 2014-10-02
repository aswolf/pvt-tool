% MieGrunDebyeHotEos - calculate Debye thermal contribution to P, E, etc
function [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,debyeDerivsFun)

    CvFac = pHotEos(end);

    % T in K
    % V in Ang^3
    % P in GPa
    R = 8.3144621/.6022*1e-3; 
    Cvm = 3*R*Natom*CvFac;

    [Tdeb0,gam0,dgamdV0] = debyeDerivsFun(V0,V0,pHotEos(1:end-1));

    % Gruneisen parameter is independent of Temp!
    [Tdeb,gam,dgamdV] = debyeDerivsFun(V,V0,pHotEos(1:end-1));
    gammaHot = gam;
    TdebyeHot= Tdeb;

    % debFunVal0 is calc at T=T0 at the *same* volume
    debFunVal  = debyeFun(Tdeb./T);
    debFunVal0 = debyeFun(Tdeb./T0);

    EHot = Cvm.*(T.*debFunVal - T0.*debFunVal0);
    CvHot = Cvm*(4*debFunVal- 3*(Tdeb./T)./(exp(Tdeb./T)-1));

    PHot = gam./V.*EHot;

    dEHotdV = 3*PHot-3*Cvm*gam.*Tdeb./V.*...
        (1./(exp(Tdeb./T)-1) - 1./(exp(Tdeb./T0)-1));

    KTHot = (gam./V-dgamdV).*EHot - gam.*dEHotdV;
end
