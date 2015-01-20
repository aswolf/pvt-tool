function [VLvl,Pval]=calcPVTBnds(Tval,PBnds,eosMod)

    bndLvl = normcdf([-1,0,1]);
    Nsamp = 30;
    Ndraw = 100;

    V0 = eosMod.pEos(1);
    evalP = @(V,pEos)(evalPressEos(pEos,eosMod,V,Tval));

    VLoP = fzero(@(V)(evalP(V,[])-PBnds(1)), V0);
    VHiP = fzero(@(V)(evalP(V,[])-PBnds(2)),.8*V0);

    VAdj = .2*abs(VLoP-VHiP);
    VBnds(1) = VHiP-VAdj;
    VBnds(2) = VLoP+VAdj;
    Vmod = linspace(VBnds(1),VBnds(2),100)';

    Pval = linspace(PBnds(1),PBnds(2),Nsamp)';

    pEosDraw = drawRandEos(eosMod,Ndraw);
    PDraw = zeros(length(Vmod),Ndraw);
    VDraw = zeros(length(Pval),Ndraw);
    for(i=1:Ndraw)
        ipEosDraw = pEosDraw(i,:);
        iPDraw = evalP(Vmod,ipEosDraw);
        PDraw(:,i) = iPDraw;
        VDraw(:,i) = interp1(iPDraw, Vmod, Pval, 'spline');
    end
    VLvl = quantile(VDraw', bndLvl);
end
