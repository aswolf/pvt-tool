function PVTeval = runPVTtool(inputFileNm)
    PVTeval = [];
    try
        [keys,sections,subsections] = inifile(inputFileNm,'readall');
    catch err
        disp('ERROR reading inputFile...');
        disp(err.message);
    end
    % NOTE: should check that keys are unique
    [data_dir,inputbase,inputext] = fileparts(inputFileNm);

    [PVTdata] = thisLoadPVTdata(keys, data_dir);
    PVTeval = thisLoadPVTeval(keys,PVTdata);

    % Perform fit
    section = 'model';
    fit_data = str2num(getKeyVal('fit_data',section,keys));
    fit_err_mod = str2num(getKeyVal('fit_err_mod',section,keys));
    if(fit_data & fit_err_mod)
        PVTeval = fitSampModPVTeval(PVTeval,[],[]);
        PVTeval = fitErrModPVTeval(PVTeval,[]);
        PVTeval = fitSampModPVTeval(PVTeval,[],[]);
    elseif(fit_data)
        PVTeval = fitSampModPVTeval(PVTeval,[],[]);
    elseif(fit_err_mod)
        disp('NOTE: this option is untested!');
        PVTeval = fitErrModPVTeval(PVTeval,[]);
    end

    % Write outputs
    section = 'output';
    update_data_table_filename = getKeyVal('update_data_table_filename',section,keys);
    eos_table_filename = getKeyVal('eos_table_filename',section,keys);
    matfilename = getKeyVal('matfilename',section,keys);

    %if(~isempty(update_data_table_filename))
    %    filepath = fullfile(data_dir,update_data_table_filename);
    %    [tblOutput,colHeader,markdownTblSep] = ...
    %        writePVTdataTbl(PVTeval.PVTdataList,filepath);
    %end
    if(~isempty(eos_table_filename))
        filepath = fullfile(data_dir,eos_table_filename);
        eosModList = [PVTeval.sampEosFit,PVTeval.sampEosPrior];
        eosModLbl = {'fit','prior'};
        [tblOutput,colHeader,markdownTblSep] = writeEosTbl(eosModList, ...
            eosModLbl,filepath,'decimal');
    end
    if(~isempty(matfilename))
        filepath = fullfile(data_dir,matfilename);
        save(filepath);
    end

    %pvt_fig_filename=
    %reduced_isotherm_fig_filename=
    %hist_fig_filename=
end
function PVTeval = thisLoadPVTeval(keys,PVTdata)
    section = 'model';
    model_name = getKeyVal('name',section,keys);
    model_material = getKeyVal('material',section,keys);
    model_fixflag = str2num(getKeyVal('fixflag',section,keys));
    model_prior_avg = str2num(getKeyVal('prior_avg',section,keys));
    model_prior_err = str2num(getKeyVal('prior_err',section,keys));
    assert(length(model_prior_avg)==length(model_prior_err),...
        'prior must define error for each value');
    assert(length(model_prior_avg)==length(model_fixflag),...
        'must define fixflag value for each parameter');
    robust_fit = str2num(getKeyVal('robust_fit',section,keys));
    robust_norm_param = str2num(getKeyVal('robust_norm_param',section,keys));
    %fitcov = str2num(getKeyVal('fitcov',section,keys));

    temp0 = str2num(getKeyVal('temp0',section,keys));
    natom = str2num(getKeyVal('natom',section,keys));
    cold_eos_fun_nm = getKeyVal('cold_eos_fun',section,keys);
    assert(strcmp(cold_eos_fun_nm,'VinetEos'),'Only VinetEos supported currently')
    cold_eos_fun = str2func(cold_eos_fun_nm);
    hot_eos_fun = str2func(getKeyVal('hot_eos_fun',section,keys));
    debye_derivs_fun = str2func(getKeyVal('debye_derivs_fun',section,keys));

    NpCold = 3;

    % Construct prior Eos
    model_prior_cov = diag(model_prior_err.^2);
    hotExtraInputs = {natom, debye_derivs_fun};
    addedThermPressFun = [];
    eosPrior = initEos('prior',model_material,model_prior_avg,model_prior_cov,...
        temp0,NpCold,cold_eos_fun,hot_eos_fun,hotExtraInputs,addedThermPressFun);

    opt = [];
    opt.robustFit =robust_fit;
    opt.robustNormParam = robust_norm_param;
    %opt.fitcov = fitcov;
    PVTeval = initPVTeval(model_name,PVTdata,eosPrior,opt);
    PVTeval.fixFlag = model_fixflag;
end
function [PVTdata] = thisLoadPVTdata(keys, data_dir)
    section = 'data';
    data_filename = getKeyVal('filename',section,keys);
    assert(~isempty(data_filename),'data filename must be provided');

    data_path = fullfile(data_dir,data_filename);
    data_name = getKeyVal('name',section,keys);
    data_grp_lbl = strsplit(getKeyVal('grp_lbl',section,keys),' ');
    err_mode = getKeyVal('err_mode',section,keys);
    mark_lbl = getKeyVal('mark_lbl',section,keys);
    mark_eos_name = getKeyVal('mark_eos',section,keys);
    try
        mark_eos = feval(['getEos_' mark_eos_name]);
    catch err
        disp('ERROR loading mark eos ...');
        disp(err.message);
    end

    % Set the options later
    opt0 = [];
    PVTdata = readPVTdataTable(data_path,data_name,data_grp_lbl,...
        err_mode,mark_lbl,mark_eos,opt0);
end
function keyVal = getKeyVal(keyname,sectionname,keys)
    isSection = strcmp(sectionname,keys(:,1));
    keyVal = keys{isSection & strcmp(keyname,keys(:,3)),4};
end
