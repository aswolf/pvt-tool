function [valMprint errMprint errSigFigMprint strFormatM] = ...
        formatTblOutput(valM,errM,numTyp,numSigFigs)
    if(isempty(numTyp))
        numTyp = repmat('f',size(valM));
    end
    if(isempty(numSigFigs))
        numSigFigs = 2;
    end

    assert(numSigFigs==2,['Custom sigfigNums is not yet implemented. '...
        'Str format output is NOT correct for '...
        'any number of sigfigs except 2'])

    powSigFig = 10.^(numSigFigs-1);
    % Round to 2 significant figures
    valSclM = 10.^floor(log10(errM));
    %errSigFigMprint = powSigFig.*roundn(errM./valSclM,1-numSigFigs);

    errMprint = valSclM.*roundn(errM./valSclM,1-numSigFigs);
    valMprint = valSclM.*roundn(valM./valSclM,1-numSigFigs);
    errSigFigMprint = powSigFig.*roundn(errM./valSclM,1-numSigFigs);

    valMprint(isnan(valMprint)) = valM(isnan(valMprint));

    % Restore zero, inf, and nan values
    errMprint(errM==0) = 0;
    errMprint(isnan(errM)) = NaN;
    errMprint(isinf(errM)) = Inf;

    errSigFigMprint(errM==0) = 0;
    errSigFigMprint(isnan(errM)) = NaN;
    errSigFigMprint(isinf(errM)) = Inf;

    % Get digit number
    %numSigFigs

    %preErrDigNum = floor(log10(valMprint))-floor(log10(errMprint));
    %totDigNum = preErrDigNum+numSigFigs;
    %postDecDigNum = totDigNum-floor(log10(valMprint))-1


    postDecDigNum = log10(valSclM)-(numSigFigs-1);
    postDecDigNum = log10(valSclM)-1;
    postDecDigNum = -min(postDecDigNum,0);
    postDecDigNum(isinf(postDecDigNum)) = 0;

    preDecDigNum = floor(log10(valM))+1;

    totDigNum = preDecDigNum+postDecDigNum+1;
    strFormatM = cell(size(valM));
    for(i=1:size(valM,1))
        for(j=1:size(valM,2))
            strFormatM{i,j} = ['%' num2str(totDigNum(i,j)) '.' ...
                num2str(postDecDigNum(i,j)) 'f'];
        end
    end
end
