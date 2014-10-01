% debyePowerLaw - return [Tdebye, gamma, dgammadV] for power law form
function [Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,peosHot)
    Tdeb0 = peosHot(1);
    gam0  = peosHot(2);
    q     = peosHot(3);

    gam    = gam0*(V./V0).^q;
    Tdeb   = Tdeb0*exp(-(gam-gam0)/q);
    dgamdV = q.*gam./V;
end
