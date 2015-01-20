function [Texp,thmExpBnds,VExpBnds]=calcThmExpBnds(Pconst,Tfoot,Tstop,dT,eosMod)
    bndLvl = normcdf([-1,0,1]);
    Ndraw = 200;

    pEosDraw = drawRandEos(eosMod,Ndraw);

    [Texp0,thmExp0,Vexp0]=calcThmExp(Pconst,Tfoot,Tstop,dT,eosMod,[]);
    thmExpDraw = zeros(Ndraw,length(Texp0));

    for(i=1:Ndraw)
        ipEos = pEosDraw(i,:);
        
        [Texp,ithmExp,iVexp,iKexp,igamexp]=calcThmExp(Pconst,Tfoot,Tstop,dT,...
            eosMod,ipEos);
        thmExpDraw(i,:) = ithmExp;
        VexpDraw(i,:) = iVexp;
    end

    thmExpBnds = quantile(thmExpDraw,bndLvl);
    VExpBnds = quantile(VexpDraw,bndLvl);
end
