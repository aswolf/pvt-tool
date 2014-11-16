function [tblOutput,colHeader,markdownTblSep] = writePVTdataTbl(PVTdata,fileNm)
    delim = '|';
    if(~exist('fileNm'))
        fileNm = '';
    end

    errMode  = PVTdata.errMode;
    measGrpID = PVTdata.measGrpID;
    uniqID = unique(measGrpID);
    IDmaxLen = 0;
    for(i=1:length(uniqID))
        ilen = length(uniqID{i});
        if(ilen > IDmaxLen)
            IDmaxLen = ilen;
        end
    end



    switch(errMode)
        case 'std'
            assert(false,'std err mode not yet implemented.');
        case 'mark'
            P  = PVTdata.Pmark;
            T  = PVTdata.T;
            V  = PVTdata.V;
            Vmark  = PVTdata.Vmark;

            TErr  = PVTdata.TErr;
            VErr  = PVTdata.VErr;
            VmarkErr  = PVTdata.VmarkErr;
            PErrTot = PVTdata.PErrTot;

            strformat = [delim '%' num2str(IDmaxLen) 's ' delim ...
                '%6.2f ' delim '%5.2f ' delim '%6.1f ' delim '%5.1f ' delim...
                '%6.2f ' delim '%5.2f ' delim '%6.2f ' delim '%5.2f' delim];
            colHeader = [delim 'grpID ' delim 'P ' delim 'PTotErr ' delim ...
                'T ' delim 'TErr ' delim 'V ' delim 'VErr ' delim ...
                'VMark ' delim 'VMarkErr' delim];
            markdownTblSep = [delim ':---' delim '---:' delim '---:' delim ...
                '---:' delim '---:' delim '---:' delim '---:' delim ...
                '---:' delim '---:' delim];
            clear tblOutput;
            for(i=1:length(V))
                istr = sprintf(strformat,measGrpID{i},P(i),PErrTot(i),...
                    T(i),TErr(i), V(i),VErr(i), Vmark(i), VmarkErr(i));
                tblOutput(i,:) = istr;
            end

    end

    if(~isempty(fileNm))
        fid = fopen(fileNm,'w');
        fprintf(fid,[colHeader '\n']);
        fprintf(fid,[markdownTblSep '\n']);
        for(i=1:size(tblOutput,1))
            fprintf(fid,[tblOutput(i,:) '\n']);
        end
        fclose(fid);
    end
end
