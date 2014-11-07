function optDefault = getPVTdataDefaultOpt();
    optDefault.calcPmark = true;
    optDefault.invVmark  = false;
    optDefault.pressFixed = false;

    optDefault.relErrPmark = false;
    optDefault.relErrVmark = false;
    optDefault.relErrV     = false;
    optDefault.relErrT     = false;

    optDefault.PTOLABS     = 0.1;

    optDefault.errModVmarkPrior = 0.3;
    optDefault.errModVPrior     = 0.3;
    optDefault.errModTPrior     = 0.3;
end
