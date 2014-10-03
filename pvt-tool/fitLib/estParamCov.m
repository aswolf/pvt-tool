% estParamCov - estimate param. covariance from numerical hess of -LogP fun
%
%  This method should finally impliment a good numerical covariance matrix solver
%    This can be broken into a number of 1D problems:
%
%    *For each parameter one-at-a-time:
%      (1) Find dp necessary for 1-sigma equiv change in cost function val
%          NOTE: target cost fun value depends on Num of dimensions
%      (2) Verify that asymmetry is within acceptable levels (~10% generally)
%    *For every pair of parameters:
%      (1) displace from min along diagonal assuming uncorrelated params
%      (2) if cost less than expected: this is the major axis direction,
%          otherwise: switch to orthogonal direction for major axis
%      (3) Find major-axis factor necessary to match target cost fun value
%      (4) Verify that asymmetry of major axis is within acceptable levels
% Write additional method to test cost multivariate normal approximation by
%  monte-carlo sampling of cov matrix and comparing Gaussian approximate
%  cost function to true cost-function.
% Using modified importance sampling, check that the credible region contains
% the right probability to within tolerance (~10%)
% If not, this implies that MCMC methods are needed
function [pfitcov,pfit] = estParamCov(pfit,priorcov,nLogPFun)
    % priorcov provides an estimate of the scale of variations
    %   - Also it shows which parameters are fixed by the prior
    psclInit = calcParamScaleInit(pfit, priorcov);
    [pfitcov,

end
% NOTE: This should use MVN slicing to determine length scale (incorporate
% correlation)
function psclInit = calcParamScaleInit(pfit, pcov)
    pSclFacDefault = 0.1;

    pscl = sqrt(diag(pcov));
    indInf = find(isinf(pscl));

    pscl(indInf) = pSclFacDefault*pfit(indInf);
end
