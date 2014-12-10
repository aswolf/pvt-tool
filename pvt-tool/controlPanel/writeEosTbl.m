function [tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList,...
        eosModLbl,fileNm,dispMode,showFlag)
    delim = '|';
    if(~exist('fileNm'))
        fileNm = '';
    end
    assert(any(strcmp(dispMode,{'paren','decimal'})), ...
        'dispMode must be either paren or decimal');
    assert(strcmp(dispMode,'decimal'),'paren mode is not yet implemented');

    pEos = eosModList(1).pEos;
    pEosM = zeros(length(pEos),length(eosModList));
    pEosCredWidM = zeros(length(pEos),length(eosModList));
    for(i=1:length(eosModList))
        ipEos = eosModList(i).pEos;
        pEosCov = eosModList(i).pEosCov;
        % if covar matrix is empty (undetermined) set values to Inf 
        %    for Null output
        if(isempty(pEosCov))
            pEosCov = diag(Inf*ones(size(ipEos)));
        end
        pEosCredWid = sqrt(diag(pEosCov));
        pEosM(:,i) = ipEos(:);
        pEosCredWidM(:,i) = pEosCredWid(:);
    end

    valM = pEosM;
    errM = pEosCredWidM;
    numTyp = [];
    numSigFigs = 2;
    [valMprint errMprint errSigFigMprint strFormatM] = ...
        formatTblOutput(valM,errM,numTyp,numSigFigs);


    paramLbl = {'$V_0$','$K_0$','$K_0''$','$\\Theta_D$','$\\gamma_0$','$q$'};

    markdownTblSep = [delim ':---' delim] ;

    switch dispMode
        case 'decimal'
            colHeader{1} = [delim delim];
            colHeader{2} = [delim delim];
            for(i=1:length(eosModList))
                colHeader{1} = [colHeader{1} eosModLbl{i} delim delim];
                colHeader{2} = [colHeader{2} 'val' delim 'err' delim];
                markdownTblSep = [markdownTblSep '---:' delim '---:' delim];
            end
            for(i=1:6)
                istr = [delim paramLbl{i} delim];
                for(j=1:length(eosModList))
                    ijVal = valMprint(i,j);
                    ijErr = errMprint(i,j);
                    ijstrFormat = [strFormatM{i,j} delim];
                    if(isinf(ijErr))
                        istr = [istr '-' delim '-' delim];
                    elseif(isnan(ijErr) | ijErr==0)
                        istr = [istr sprintf(ijstrFormat,ijVal) '-' delim];
                    else
                        istr = [istr sprintf(ijstrFormat,ijVal) ...
                            sprintf(ijstrFormat,ijErr) ];
                    end
                end
                tblOutput{i} = istr;
            end
        case 'paren'
            colHeader{1} = [delim delim];
            for(i=1:length(eosModList))
                colHeader{1} = [colHeader{1} eosModLbl{i} delim];
                markdownTblSep = [markdownTblSep '---:' delim];
            end
            for(i=1:6)
                istr = [delim paramLbl{i} delim];
                for(j=1:length(eosModList))
                    ijVal = valMprint(i,j);
                    ijErrFig = errSigFigMprint(i,j);
                    ijstrFormat = strFormatM{i,j};
                    if(isinf(ijErrFig))
                        istr = [istr '-' delim];
                    elseif(isnan(ijErrFig) | ijErrFig==0)
                        istr = [istr sprintf(ijstrFormat,ijVal) delim];
                    else
                        istr = [istr sprintf(ijstrFormat,ijVal) ...
                            '(' sprintf('%d',ijErrFig) ')' delim];
                    end
                end
                tblOutput{i} = istr;
            end
    end

    if(~isempty(fileNm))
        fid = fopen(fileNm,'w');
        fprintf(fid,[colHeader{1} '\n']);
        fprintf(fid,[markdownTblSep '\n']);
        if(strcmp(dispMode,'decimal'))
            fprintf(fid,[colHeader{2} '\n']);
        end
        for(i=1:length(tblOutput))
            fprintf(fid,[tblOutput{i} '\n']);
        end
        fclose(fid);
    end
end
