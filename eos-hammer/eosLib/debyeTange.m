% debyeTange - return [Tdebye, gamma, dgammadV] using form from Tange 2009
function [Tdeb, gam, dgamdV] = debyeTange(V,V0,peosHot)
    Tdeb0 = peosHot(1);
    gam0  = peosHot(2);
    a     = peosHot(3);
    b     = peosHot(4);

    gam    = gam0*(1+ a*((V./V0).^b-1));
    Tdeb   = Tdeb0*(V./V0).^(-(1-a)*gam0).*exp(-(gam-gam0)/b);
    dgamdV = gam0/V0*a*b*(V./V0).^(b-1);
end
