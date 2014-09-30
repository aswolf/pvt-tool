% debyeFun - calculates debye function for vector of values using for loop
%
% cannot handle values near realmin... should never happen for EOS applications
function Dvec = debyeFun(x)
    Dvec = ones(size(x));
    for(i=1:length(x))
        xi = abs(x(i));
        if xi >= realmin
            D1 = 3*quadgk(@debye1_integrand,0,xi)./xi.^3;
        else
            D1 = 1;
        end
        Dvec(i) = D1;
    end
end
function y = debye1_integrand(t)
    y = ones(size(t));
    nz = (abs(t) >= realmin);
    y(nz) = t(nz).^3./(exp(t(nz))-1);
end
