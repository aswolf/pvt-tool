% debyePowerLaw - return [Tdebye, gamma, dgammadV] for power law form
function [Tdeb, gam, dgamdV] = debyeShiftPowerLaw(V,V0,peosHot)
    Tdeb0 = peosHot(1);
    gam1  = peosHot(2);
    q     = peosHot(3);

    relV = V./V0;

    gam    = gam1*relV.^q + 0.5;
    Tdeb   = Tdeb0*relV.^(-0.5).*exp(gam1*(1-relV.^q)/q);
    dgamdV = q*gam1/V0*relV.^(q-1);
end
