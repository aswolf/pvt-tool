% MieGrunDebyeHotEos - calculate Debye thermal contribution to P, E, etc
function [Phot,Ehot,Cvhot,KThot,Tdebyehot,gammahot] = ...
        MieGrunDebyeHotEos(V,T,V0,T0,Natom,peosHot,debyeTderivsFun)

    CvFac = peosHot(end);

    % T in K
    % V in Ang^3
    % P in GPa
    R = 8.3144621/.6022*1e-3; 
    Cvm = 3*R*Natom*CvFac;

    [Tdeb0,gam0,dgamdV0] = debyeTderivsFun(V0,V0,peosHot(1:end-1));

    % Gruneisen parameter is independent of Temp!
    [Tdeb,gam,dgamdV] = debyeTderivsFun(V,V0,peosHot(1:end-1));
    gammahot = gam;
    Tdebyehot= Tdeb;

    % debFunVal0 is calc at T=T0 at the *same* volume
    debFunVal  = debyeFun(Tdeb./T);
    debFunVal0 = debyeFun(Tdeb./T0);

    Ehot = Cvm.*(T.*debFunVal - T0.*debFunVal0);
    Cvhot = Cvm*(4*debFunVal- 3*(Tdeb./T)./(exp(Tdeb./T)-1));

    Phot = gam./V.*Ehot;

    dEhotdV = 3*Phot-3*Cvm*gam.*Tdeb./V.*...
        (1./(exp(Tdeb./T)-1) - 1./(exp(Tdeb./T0)-1));

    KThot = (gam./V-dgamdV).*Ehot - gam.*dEhotdV;
    %KThot = -Cvm*((dgamdV-gam./V).*(T.*debFunVal-T0*debFunVal0) ...
    %    + 3*gam.^2./V.*(T.*debFunVal - Tdeb./(exp(Tdeb./T)-1)));
end
