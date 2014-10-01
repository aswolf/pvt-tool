function [PCold,KCold,dFCold,KPCold] = VinetEos(V,pColdEos)
    assert(length(pColdEos) >= 3,'pcoldEos must have at least 3 elements');
    V0  = pColdEos(1);
    K0  = pColdEos(2);
    K0P = pColdEos(3);

    nu = 3/2*(K0P -1);
    x = (V./V0).^(1/3);
    PCold = (3*K0.*(1-x)./x.^2).*exp(nu*(1-x));

    % Check if need to calculate Bulk Modulus
    if(nargout() >= 2)
        KCold = -K0*(x.*(-nu*(1-x)+1)-2)./x.^2.*exp(nu*(1-x));
    else
        KCold = [];
    end

    % Check whether calculating helmholtz energy is necessary 
    % NOTE: units of energy are [press units]*[Vol units]
    if(nargout() >= 3)
        dFCold = 9*K0*V0/nu^2 * (1 - exp(nu*(1-x)).*(nu*(x-1)+1));
    else
        dFCold = [];
    end

    % Check if need to calculate 1st press deriv of bulk modulus
    if(nargout() >= 4)
        KPCold = -K0./(3*x.^2.*KCold).*exp(nu*(1-x)).*...
            (x.*(nu^2*x.*(x-1) + nu*(x-3)+1)-4);
    else
        KPCold = [];
    end
end
