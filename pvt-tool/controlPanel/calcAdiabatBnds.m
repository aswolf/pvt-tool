function [Pad,VadBnds,KadBnds]=calcAdiabatBnds(Tfoot,Pfoot,Pstop,dlogV,eosMod)
    bndLvl = normcdf([-1,0,1]);
    Ndraw = 200;
    Nsamp = 30;

    pEosDraw = drawRandEos(eosMod,Ndraw);

    [Pad0,Vad0,Tad0,Kad0,gamad0,thmExpad0]=...
        calcAdiabat(Tfoot,Pfoot,Pstop,dlogV,eosMod,[]);
    
    Pad  = linspace(Pstop,Pfoot,Nsamp);
    for(i=1:Ndraw)
        ipEos = pEosDraw(i,:);
        
        [iPad,iVad,iTad,iKad,igamad,ithmExpad]=...
            calcAdiabat(Tfoot,Pfoot,Pstop,dlogV,eosMod,ipEos);
        VadDraw(i,:) = interp1(iPad,iVad,Pad,'spline');
        KadDraw(i,:) = interp1(iPad,iKad,Pad,'spline');
    end

    VadBnds = quantile(VadDraw,bndLvl);
    KadBnds = quantile(KadDraw,bndLvl);
end
