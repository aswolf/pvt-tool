function [costDraw,costDrawApp,pdraw,costDrawMin] = ...
        testCostFun(Ndraw,pavg,pcov,fixFlag,costFun)

    costDrawMin = costFun(pavg);

    [pFree,pFree,pcovFree] = getFreeParams(pavg,pavg,pcov,fixFlag);

    pHessFree = inv(pcovFree);

    pdraw = zeros(Ndraw,length(pavg));
    costDraw = zeros(Ndraw,1);
    costDrawApp = zeros(Ndraw,1);

    for(i=1:Ndraw)
        ipdrawFree = mvnrnd(pFree,pcovFree,1);
        icostDrawApp = costDrawMin +...
            0.5*(ipdrawFree-pFree)*pHessFree*(ipdrawFree-pFree)';

        ipdraw = getAllParams(ipdrawFree,pavg,fixFlag);
        icostDraw = costFun(ipdraw);

        pdraw(i,:) = ipdraw;
        costDraw(i) = icostDraw;
        costDrawApp(i) = icostDrawApp;
    end
end
