function eos = initEos(T0,pColdEos,pHotEos,fullEosModTyp,coldEosTyp,hotEosTyp)

    if(nargin==0)
        eos = struct(...
        'fullMieGrunDebTyp',@fullMieGrunDebTyp,...
        'fullThermExpTyp',@fullThermExpTyp,...
        'coldVinetTyp',@coldVinetTyp,...
        'coldBirchMurnTyp',@coldBirchMurnTyp,...
        'hotDebyePowerLawTyp',@hotDebyePowerLawTyp,...
        'hotDebyeTangeTyp',@hotDebyeTangeTyp);
        return eos;
    end
    %checkInput(measGrpID,P,PErr,VMark,VMarkErr,VSamp,VSampErr,T,Terr);

    %[measGrpID,P,PErr,VMark,VMarkErr,VSamp,VSampErr,T,Terr] = ...
    %    cleanInput(measGrpID,P,PErr,VMark,VMarkErr,VSamp,VSampErr,T,Terr);


    % Initialize error model for PVT data using external method in fitLib
    initCalibState = initErrModPvt();

    % If not defined by user, all data put into 1 meas grp
    if(isempty(measGrpID))
        measGrpID = zeros(Ndat,1);
    end
    measTyp    = detectMeasTyp(P,T);


    % Initialize pvt struct
    pvt = struct(...
        'calibState', initCalibState ,...
        'measTyp'   , measTyp ,...
        'measGrpID' , measGrpID, ...
        'P'         , P ,...
        'PErr'      , PErr ,...
        'VMark'     , VMark ,...
        'VMarkErr'  , VMarkErr ,...
        'VSamp'     , VSamp ,...
        'VSampErr'  , VSampErr ,...
        'T'         , T ,...
        'TErr'      , TErr ,...
        'coldCompressTyp',@coldCompressTyp,...
        'hotCompressTyp' ,@hotCompressTyp,...
        'thermExpTyp'    ,@thermExpTyp);
end
function measTyp = detectMeasTyp(P,T)
    pvt = initPvtData();
    %    (1) cold compression data
    %    (2) ambient press therm exp data
    %    (3) general high P-T data
    T0 = 300;
    P0 = 0;

    TOL_P = 1;
    TOL_T = 10

    isColdCompress = (abs(T-T0) < TOL_T);
    isThermExp     = (abs(P-P0) < TOL_P & ~isColdCompress);

    measTyp = ones(length(P),1);

    % cold compression data is in group 1
    measTyp(isColdCompress) = pvt.coldCompressTyp;
    % thermal expansion data is in group 2
    measTyp(isThermExp) = pvt.thermExpTyp;
    % General high P-T data is in group 3
    measTyp(~isColdCompress & ~isThermExp) = pvt.hotCompressTyp;
end

function [measGrpID,P,PErr,VMark,VMarkErr,VSamp,VSampErr,T,Terr] = ...
        cleanInput(measGrpID,P,PErr,VMark,VMarkErr,VSamp,VSampErr,T,Terr)

    Ndat = length(VSamp);

    measGrpID = measGrpID(:);
    P        = P(:);
    PErr     = PErr(:);
    VMark    = VMark(:);
    VMarkErr = VMarkErr(:);
    VSamp    = VSamp(:);
    VSampErr = VSampErr(:);
    T        = T(:);
    Terr     = Terr(:);

    makeFullArr = @(arr)(arr*ones(Ndat,1));

    if(length(measGrpID)==1) measGrpID = makeFullArr(measGrpID); end
    if(length(P)==1)           P           = makeFullArr(P); end
    if(length(PErr)==1)        PErr        = makeFullArr(PErr); end
    if(length(VMark)==1)       VMark       = makeFullArr(VMark); end
    if(length(VMarkErr)==1)    VMarkErr    = makeFullArr(VMarkErr); end
    if(length(VSamp)==1)       VSamp       = makeFullArr(VSamp); end
    if(length(VSampErr)==1)    VSampErr    = makeFullArr(VSampErr); end
    if(length(T)==1)           T           = makeFullArr(T); end
    if(length(TErr)==1)        TErr        = makeFullArr(TErr); end
end
function checkInput(measGrpID,P,PErr,VMark,VMarkErr,VSamp,VSampErr,T,Terr)
    assert(~isempty(VSamp),'VSamp cannot be empty.');
    Ndat = length(VSamp);

    % Verify that we have enough data to determine P,VSamp,T for every datapt
    assert(length(T) == 1 | length(T) == Ndat, ...
        'T must be either a constant or provided for every datapt.');
    if(isempty(VMark))
        assert(length(P) == Ndat, ...
            'If VMark not provided, P must be provided for every datapt.');
    else
        assert(length(VMark) == Ndat, ...
            'If provided, VMark must be given for every datapt.');
    end

    assert(length(measGrpID) == 1 | length(measGrpID) == Ndat, ...
        'measGrpID must be either a constant or provided for every datapt.');
    assert(isempty(PErr) | length(PErr) == 1 | length(PErr) == Ndat, ...
        'PErr must be empty, a constant, or provided for every datapt.');
    assert(isempty(TErr) | length(TErr) == 1 | length(TErr) == Ndat, ...
        'TErr must be empty, a constant, or provided for every datapt.');
    assert(isempty(VSampErr) | length(VSampErr) == 1 | length(VSampErr) == Ndat, ...
        'VSampErr must be empty, a constant, or provided for every datapt.');
    assert(isempty(VMarkErr) | length(VMarkErr) == 1 | length(VMarkErr) == Ndat, ...
        'VMarkErr must be empty, a constant, or provided for every datapt.');
end
function typID = coldCompressTyp()
    typID = 0;
end
function typID = hotCompressTyp()
    typID = 1;
end
function typID = thermExpTyp()
    typID = 2;
end
