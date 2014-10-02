function [PVTmodstruc] = evalEosPVTwrap(V,T,peos,T0,eosColdFun,eosMGDFun)

    peos    = eosStruc.eos;
    T0      = eosStruc.T0;
    eosColdFun = eosStruc.eosColdFun;
    eosMGDFun  = eosStruc.eosMGDFun;

    eosHotFun = @(V,T,peos)(eosColdFun(V,peos(1:3))+...
        eosMGDFun(V,T,peos(1),T0,peos(4:end)));

    [Pmod,dPdVmod,dPdTmod   ] = evalEosPVT(V,    T,peos,eosHotFun);
    [Pmark,dPdVmark,dPdTmark] = evalEosPVT(Vmark,T,peos,eosHotFun);

