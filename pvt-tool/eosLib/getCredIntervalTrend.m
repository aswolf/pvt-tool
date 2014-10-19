function [trendBnds,trendMed,trend] = getCredIntervalTrend(pEos,pcovEos,...
        trendFun,Ndraw,pLvl)
    assert(pLvl<0.5 & pLvl > 0,'pLvl must be less than .5 and greater than 0');

    trend = trendFun(pEos);
    % Keep fixed parameters, marked by Nans, as fixed
    pcovEos(isnan(pcovEos)) = 0;
    pdrawEos = mvnrnd(pEos,pcovEos,Ndraw);

    for(i=1:Ndraw)
        itrend = trendFun(pdrawEos(i,:));
        if(i==1)
            trendDraw = zeros(Ndraw,length(itrend));
        end
        trendDraw(i,:) = itrend;
    end
    trendBnds = quantile(trendDraw,.5+pLvl*[-1 1]);
    trendMed = quantile(trendDraw,.5);
end



