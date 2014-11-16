% formatValErrStr - format Value/Error pair into string for use in Tables
function [valStr errStr errSigFigStr] = formatValErrStr(val,err,numSigFigs,numTyp)
    if(isempty(numSigFigs))
        numSigFigs = 2;
    end
    if(isempty(numTyp))
        numTyp = 'float';
    end

    orderVal = floor(log10(abs(val)));
    orderErr = floor(log10(abs(err)));

    orderDiff = orderVal-orderErr;
    orderTotRange = orderDiff+numSigFigs;

    if(strcmp(numTyp,'float'))
        % set field width to total order range plus 1 for decimal character
        %  unless total order range is less than orderVal
        if(orderTotRange > orderVal+1)
            fieldWidVal = orderTotRange + 1;
            decimalPrecVal = orderTotRange - (orderVal+1);
        else
            fieldWidVal = orderVal + 1;
            decimalPrecVal = 0;
        end

        % Error formating for diff cases
        if(orderErr < 0)
            % error less than 1
            fieldWidErr =  abs(orderErr) + (numSigFigs) + 1;
            decimalPrecErr = fieldWidErr - 2;
        else
            if(orderErr+1-numSigFigs >= 0)
                % last significant fig is still greater than fraction
                fieldWidErr = orderErr + 1;
                decimalPrecErr = 0;
            else
                % error spans decimal point 
                fieldWidErr = numSigFigs+1;
                decimalPrecErr = -(orderErr+1-numSigFigs);
            end
        end
    elseif(strcmp(numTyp,'int'))
        fieldWidVal = orderVal + 1;
        decimalPrecVal = 0;
        fieldWidErr = orderErr + 1;
        decimalPrecErr = 0;
    elseif(strcmp(numTyp,'exp'))
        fieldWidVal = orderTotRange + 1;
        decimalPrecVal = orderTotRange-1;
        fieldWidErr = numSigFigs + 1;
        decimalPrecErr = numSigFigs - 1;
    end
    % If negative, add to field width to allow for negative sign
    if(val<0)
        fieldWidVal = fieldWidVal+1;
    end

    strFormatVal = ['%' num2str(fieldWidVal) '.' num2str(decimalPrecVal)];
    strFormatErr = ['%' num2str(fieldWidErr) '.' num2str(decimalPrecErr)];
    strFormatErrSigFig = ['%' num2str(numSigFigs) '.f'];

    if(any(strcmp(numTyp,{'float','int'})))
        strFormatVal = [strFormatVal 'f'];
        strFormatErr = [strFormatErr 'f'];
    elseif(strcmp(numTyp,'exp'))
        strFormatVal = [strFormatVal 'e'];
        strFormatErr = [strFormatErr 'e'];
    end
    valStr = sprintf(strFormatVal,val);
    errStr = sprintf(strFormatErr,err);
    errSigFigStr = sprintf(strFormatErrSigFig,err/10.^(orderErr-numSigFigs+1));
end
