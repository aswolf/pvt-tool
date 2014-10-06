% getEstParamDefaultOpt - get default options for estParamLib
function optDefault = getEstParamDefaultOpt()
    % Assume standard least squares fitting
    optDefault.robustFit   = false;
    % Assume standard least squares penalty for prior
    optDefault.robustPrior = false;
    % if robust fitting desired, use Student's t distribution with 
    %  reasonable normality parameter value (nu=5)
    optDefault.robustNormParam = 5;
    % use standard fminunc function for finding minimum
    optDefault.minAlgFun = @fminunc;
    optDefault.LargeScale = 'off';
    optDefault.Display = 'off';
    % repeat min operation to be sure at min
    optDefault.NfitIter  = 2;
end
