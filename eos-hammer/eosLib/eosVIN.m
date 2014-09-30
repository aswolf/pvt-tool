function [P,dE,KLoc,KPLoc] = eosVIN(V,params)
    assertTrue(length(params) >= 3,'params must have at least 3 elements');
    V0  = params(1);
    K0  = params(2);
    K0P = params(3);

    nu = 3/2*(K0P -1);
    x = (V./V0).^(1/3);
    P = (3*K0.*(1-x)./x.^2).*exp(nu*(1-x));

    % Check whether calculating energy is necessary 
    % NOTE: units of energy are [press units]*[Vol units]
    if(nargout() >= 2)
        dE = 9*K0*V0/nu^2 * (1 - exp(nu*(1-x)).*(nu*(x-1)+1));
    else
        dE = [];
    end

    % Check if need to calculate Bulk Modulus
    if(nargout() >= 3)
        KLoc = -K0*(x.*(-nu*(1-x)+1)-2)./x.^2.*exp(nu*(1-x));
    else
        KLoc = [];
    end

    % Check if need to calculate 1st press deriv of bulk modulus
    if(nargout() >= 4)
        KPLoc = -K0./(3*x.^2.*KLoc).*exp(nu*(1-x)).*...
            (x.*(nu^2*x.*(x-1) + nu*(x-3)+1)-4);
    else
        KPLoc = [];
    end
end
