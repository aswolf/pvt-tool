% estBasinHess - estimate hessian of scalar fun within depth region of minimum 
function [pHess,pAsym] = estBasinHess(basinFun, pmin, opt)
    % Whole Function to be implemented later
    %opt.pscl
    %opt.depth
    %opt.TOL
    %opt.maxDev
    %optDefault = getEstParamDefaultOpt();
%  This method should finally impliment a good numerical covariance matrix solver
%    This can be broken into a number of 1D problems:
%
%    *For each parameter one-at-a-time:
%      (1) Find dp necessary for depth change in basin function val
%          NOTE: target basin fun value is input 
%      (2) Calc asymmetry of individual param slices
%    *For every pair of parameters:
%      (1) displace from min along diagonal assuming uncorrelated params
%      (2) Fit 2D quadratic to pure samples (along x1 and x2) and mixed sample
%      (3) Jump to new point where expected depth = target depth along major
%          axis
%      (4) Iterate until convergance
%      (5) Calc asymmetry of major axis param slice by jumping to opposite
%          side of minimum
% Write additional method to test cost multivariate normal approximation by
%  monte-carlo sampling of cov matrix and comparing Gaussian approximate
%  cost function to true cost-function.
% Using modified importance sampling, check that the credible region contains
% the right probability to within tolerance (~10%)
% If not, this implies that MCMC methods are needed
    pHess = [];
    pAsym = [];

end

