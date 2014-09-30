function [Pthm,KTthm,Cvthm,gamma,delEthm] = calcPthmMGDpowlaw(V,T,V0,T0,peosHot)
    % T in K
    % V in Ang^3
    % P in GPa
    R = 8.3144621/.6022*1e-3; 

    thetaD0     = peosHot(1);
    gamma0      = peosHot(2);
    q           = peosHot(3);
    relCvmperCell = peosHot(4);


    Cvm = 3*R*relCvmperCell;

    gamma  = gamma0.*(V./V0).^q;
    thetaD = thetaD0*exp((gamma0-gamma)./q);

    debyeFunVal  = debyeFun(thetaD./T);
    debyeFunVal0 = debyeFun(thetaD0./T0);

    Pthm = Cvm*gamma./V.*(T.*debyeFunVal - T0*debyeFunVal0);

    KTthmFun=@(T)(-Cvm.*T.*gamma./V.*(debyeFunVal.*(q-1+3*gamma)...
        - 3*gamma.*(thetaD./T)./(exp(thetaD./T)-1)));
    KTthm = KTthmFun(T) - KTthmFun(T0);

    Cvthm = Cvm*(4*debyeValFun- 3*thetaD./T./(exp(thetaD./T)-1));

    delEthm = + Cvm*(T./3.*(debyeFunVal ...
        -3.0*log((1-exp(-thetaD./T))./(1-exp(-thetaD0/T0))) ...
        -4*debyeFunVal0) + T.*debyeFunVal);
end


