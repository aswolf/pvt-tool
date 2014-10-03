function tests = eosLibTest
    tests = functiontests(localfunctions);
end
function testDebyeFunVals(testCase)
    TOL = 1e-5;
    x = [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.6 1.8 ...
        2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.5  ...
        6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0];
    DxTbl = [1.000000 0.963000 0.926999 0.891995 0.857985 0.824963 0.792924 ...
        0.761859 0.731759 0.702615 0.674416 0.647148 0.620798 0.595351 ...
        0.570793 0.524275 0.481103 0.441129 0.404194 0.370137 0.338793 ...
        0.309995 0.283580 0.259385 0.237252 0.217030 0.198571 0.181737 ...
        0.166396 0.152424 0.139704 0.128129 0.117597 0.095241 0.077581 ...
        0.063604 0.052506 0.043655 0.036560 0.030840 0.026200 0.022411 ...
        0.019296];
    Dxcalc = debyeFun(x);
    calcErr = Dxcalc - DxTbl;
    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testVinetEos_P0(testCase)
    V0 = 10;
    pEos = [V0,100,4.5];
    verifyEqual(testCase,VinetEos(V0,pEos),0);
end
function testVinetEos_K0(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 100;
    pEos = [V0,K0,4.5];
    dVstep = .001*[-1,1];
    V = V0*(1+dVstep);
    P = VinetEos(V,pEos);
    K0num = -V0*diff(P)/diff(V);

    calcErr = K0num/K0-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testVinetEos_KP0(testCase)
    TOL = 1e-5;

    V0 = 10;
    K0 = 100;
    KP0 = 4.5;
    pEos = [V0,K0,KP0];
    dVstep = 1e-4*[-1:1:1];
    V = V0*(1+dVstep);
    P = VinetEos(V,pEos);

    % Perform 2nd order Taylor expansion to determine KP0
    K0num = -V0*diff(P)/diff(V);
    designM = [.5*P(:).^2 P(:) ones(size(P(:)))];
    polyv = designM\V(:);
    KP0num = V0*polyv(1)/polyv(2)^2 -1;

    calcErr = KP0num/KP0-1;
    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testVinetEos_Peqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 200;
    pEos = [V0,K0,4.5];

    V = V0*.7*(1+1e-4*[-2:2]);
    [P K dF] = VinetEos(V,pEos);

    Pnum = -central_diff(dF,V);
    calcErr = Pnum(3) - P(3);

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testVinetEos_Keqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 200;
    pEos = [V0,K0,4.5];

    V = V0*.7*(1+1e-4*[-2:2]);
    [P K] = VinetEos(V,pEos);

    Knum = -V.*central_diff(P,V);
    calcErr = Knum(3)./K(3)-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end

function testVinetEos_KPeqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 200;
    pEos = [V0,K0,4.5];

    V = V0*.7*(1+1e-4*[-2:2]);
    [P K dF KP] = VinetEos(V,pEos);

    KPnum = central_diff(K,P);
    calcErr = KPnum(3)./KP(3)-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end

%%%%%%%%%%%%%%%%%%%%%%
%  Test Power Law fun
%%%%%%%%%%%%%%%%%%%%%%
function testDebyePowerLaw_Tdeb0(testCase)
    V0 = 10;

    Tdeb0 = 500;
    gam0 = 1.7;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q];

    V = V0;
    [Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    verifyEqual(testCase,Tdeb,Tdeb0)
end
function testDebyePowerLaw_gam0(testCase)
    V0 = 10;

    Tdeb0 = 500;
    gam0 = 1.7;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q];

    V = V0;
    [Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    verifyEqual(testCase,gam,gam0)
end
function testDebyePowerLaw_TdebEqn(testCase)
    TOL = 1e-5;
    V0 = 10;

    Tdeb0 = 500;
    gam0 = 1.7;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q];

    V = V0*.7*(1+1e-4*[-2:2]);
    [Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    gamnum = -V./Tdeb.*central_diff(Tdeb,V);

    calcErr = gamnum(3) - gam(3);
    verifyTrue(testCase,abs(calcErr)<TOL)
end
function testDebyePowerLaw_dgamdVEqn(testCase)
    TOL = 1e-5;
    V0 = 10;

    Tdeb0 = 500;
    gam0 = 1.7;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q];

    V = V0*.7*(1+1e-4*[-2:2]);
    [Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);

    dgamdVnum = central_diff(gam,V);

    calcErr = dgamdVnum(3) - dgamdV(3);
    verifyTrue(testCase,abs(calcErr)<TOL)
end

%%%%%%%%%%%%%%%%%%%%%%
%  Test Tange 09 fun
%%%%%%%%%%%%%%%%%%%%%%
function testDebyeTange_Tdeb0(testCase)
    V0 = 74.698;

    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;
    pHotEos = [Tdeb0 gam0 a b];

    V = V0;
    [Tdeb, gam, dgamdV] = debyeTange(V,V0,pHotEos);
    verifyEqual(testCase,Tdeb,Tdeb0)
end
function testDebyeTange_gam0(testCase)
    V0 = 74.698;

    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;
    pHotEos = [Tdeb0 gam0 a b];

    V = V0;
    [Tdeb, gam, dgamdV] = debyeTange(V,V0,pHotEos);
    verifyEqual(testCase,gam,gam0)
end
function testDebyeTange_TdebEqn(testCase)
    TOL = 1e-5;
    V0 = 74.698;

    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;
    pHotEos = [Tdeb0 gam0 a b];

    V = V0*.7*(1+1e-4*[-2:2]);
    [Tdeb, gam, dgamdV] = debyeTange(V,V0,pHotEos);
    gamnum = -V./Tdeb.*central_diff(Tdeb,V);

    calcErr = gamnum(3) - gam(3);
    verifyTrue(testCase,abs(calcErr)<TOL)
end
function testDebyeTange_dgamdVEqn(testCase)
    TOL = 1e-5;
    V0 = 74.698;

    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;
    pHotEos = [Tdeb0 gam0 a b];

    V = V0*.7*(1+1e-4*[-2:2]);
    [Tdeb, gam, dgamdV] = debyeTange(V,V0,pHotEos);

    dgamdVnum = central_diff(gam,V);

    calcErr = dgamdVnum(3) - dgamdV(3);
    verifyTrue(testCase,abs(calcErr)<TOL)
end


%%%%%%%%%%%%%%%%%%%%%%
%  Test MieGrunDebye Hot Eos fun
%%%%%%%%%%%%%%%%%%%%%%
function testMieGrunDebyeHotEos_PHot0(testCase)
    V0 = 10;
    T0 = 300;
    Natom = 4;
    Tdeb0 = 500;
    gam0 = 1.5;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q 1.0];

    V = V0;
    T = T0;
    %[Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,@debyePowerLaw);

    verifyEqual(testCase,PHot,0);
end
function testMieGrunDebyeHotEos_EHot0(testCase)
    V0 = 10;
    T0 = 300;
    Natom = 4;
    Tdeb0 = 500;
    gam0 = 1.5;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q 1.0];

    V = V0;
    T = T0;
    %[Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,@debyePowerLaw);

    verifyEqual(testCase,EHot,0);
end
function testMieGrunDebyeHotEos_CvEqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    T0 = 300;
    Natom = 4;
    Tdeb0 = 500;
    gam0 = 1.5;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q 1.0];

    V = V0*.7;
    T = Tdeb0*(1+1e-4*[-2:2]);

    %[Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,@debyePowerLaw);

    Cvnum = central_diff(EHot,T);
    calcErr = Cvnum(3)./CvHot(3)-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testMieGrunDebyeHotEos_KHotEqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    T0 = 300;
    Natom = 4;
    Tdeb0 = 500;
    gam0 = 1.5;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q 1.0];

    V = V0*.7*(1+1e-4*[-2:2]);
    T = Tdeb0;

    %[Tdeb, gam, dgamdV] = debyePowerLaw(V,V0,pHotEos);
    [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,@debyePowerLaw);

    KHotnum = -V.*central_diff(PHot,V);
    calcErr = KHotnum(3)./KTHot(3)-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end

%%%%%%%%%%%%%%%%%%%%%%
%  Test MieGrunEinstein Hot Eos fun
%%%%%%%%%%%%%%%%%%%%%%
function testMieGrunEinsteinHotEos_PHot0(testCase)
    TOL = 1e-5;
    V0 = 10;
    T0 = 300;
    Natom = 4;
    Tdeb0 = 500;
    gam0 = 1.5;
    q = 1.2;
    pHotEos = [Tdeb0 gam0 q 1.0];

    V = V0*.7*(1+1e-4*[-2:2]);
    T = Tdeb0;
    [PHot,KTHot,CvHot,gammaHot,EHot,TdebyeHot] = ...
        MieGrunEinsteinHotEos(V,T,V0,T0,Natom,pHotEos,@debyePowerLaw);
    verifyTrue(testCase,false);
end

%%%%%%%%%%%%%%%%%%%%%%
%  Test ThermAddEos (thermal Pressure form)
%%%%%%%%%%%%%%%%%%%%%%
function testCalcPressThermAddEos_Tange2012(testCase)
    T0 = 300;
    V0 = 162.373;
    K0 = 258.4;
    KP0= 4.10;

    Natom = 4*5;
    Tdeb0 = 940;
    gam0  = 1.55;
    q = 1.1;
    
    Tcol = [300 1000 2000 3000 4000];
    Vrow = V0*[1.00 0.98 0.96 0.94 0.92 0.90 0.88 0.86 0.84 0.82 0.80 0.78 0.76];
    [T,V] = meshgrid(Tcol,Vrow);


    TOL = 2e-2;
    PTbl = [...
          0.00    4.83   12.57   20.42   28.30;...
          5.44   10.22   17.93   25.77   33.63;...
         11.46   16.20   23.88   31.69   39.53;...
         18.13   22.82   30.47   38.26   46.09;...
         25.52   30.15   37.77   45.54   53.35;...
         33.70   38.27   45.86   53.61   61.40;...
         42.75   47.27   54.84   62.56   70.33;...
         52.79   57.25   64.78   72.49   80.23;...
         63.92   68.32   75.82   83.50   91.22;...
         76.28   80.61   88.08   95.73  103.43;...
         89.99   94.26  101.69  109.32  117.00;...
        105.23  109.44  116.83  124.43  132.09;...
        122.19  126.32  133.68  141.25  148.89;...
        ];

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 q 1.0];
    
    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyePowerLaw;
    hotEosFun  = @(V,T,V0,T0,pHotEos)...
        (MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,debyeDerivsFun));
    elecThermPressFun = [];

    [P,KT,Cv,gam] = calcPressThermAddEos(V(:),T(:),T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,elecThermPressFun);
    PmodTbl = reshape(P,length(Vrow),[]);

    calcErr = PmodTbl-PTbl;

    verifyTrue(testCase,all(reshape(abs(calcErr),[],1) < TOL));
end

function testCalcPressThermAddEos_Tange2009(testCase)
    T0 = 300;
    V0 = 74.698;
    K0 = 160.63;
    KP0= 4.367;

    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyeTange;
    Natom = 4*2;
    Tdeb0 = 761;
    gam0  = 1.442;
    a = 0.138;
    b = 5.4;
    Tcol = [300    500   1000   1500   2000   2500   3000   4000];

    Vrow = V0*[1.15 1.14 1.13 1.12 1.11 1.10 1.09 1.08 1.07 1.06 1.05 ...
        1.04 1.03 1.02 1.01 1.00 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 ...
        0.91 0.90 0.89 0.88 0.87 0.86 0.85 0.84 0.83 0.82 0.81 0.80 0.79 ...
        0.78 0.77 0.76 0.75 0.74 0.73 0.72 0.71 0.70 0.69 0.68 0.67 0.66 ...
        0.65];

    [T,V] = meshgrid(Tcol,Vrow);

    TOL = 2e-2;
    PTbl = [...
           NaN     NaN     NaN     NaN     NaN     NaN    0.55    6.97 ; ...
           NaN     NaN     NaN     NaN     NaN     NaN    1.21    7.61 ; ...
           NaN     NaN     NaN     NaN     NaN     NaN    1.92    8.30 ; ...
           NaN     NaN     NaN     NaN     NaN     NaN    2.68    9.05 ; ...
           NaN     NaN     NaN     NaN     NaN    0.32    3.50    9.85 ; ...
           NaN     NaN     NaN     NaN     NaN    1.20    4.37   10.71 ; ...
           NaN     NaN     NaN     NaN     NaN    2.13    5.30   11.63 ; ...
           NaN     NaN     NaN     NaN     NaN    3.12    6.29   12.62 ; ...
           NaN     NaN     NaN     NaN    1.02    4.18    7.34   13.67 ; ...
           NaN     NaN     NaN     NaN    2.14    5.30    8.46   14.79 ; ...
           NaN     NaN     NaN    0.18    3.33    6.48    9.65   15.98 ; ...
           NaN     NaN     NaN    1.44    4.59    7.74   10.91   17.25 ; ...
           NaN     NaN     NaN    2.77    5.92    9.08   12.25   18.59 ; ...
           NaN     NaN    1.04    4.17    7.33   10.49   13.67   20.02 ; ...
           NaN     NaN    2.52    5.65    8.82   11.99   15.17   21.53 ; ...
          0.00    1.06    4.09    7.22   10.39   13.57   16.76   23.14 ; ...
          1.65    2.71    5.74    8.88   12.06   15.25   18.44   24.84 ; ...
          3.39    4.45    7.48   10.63   13.82   17.01   20.22   26.64 ; ...
          5.23    6.29    9.32   12.48   15.68   18.88   22.10   28.54 ; ...
          7.17    8.22   11.26   14.43   17.64   20.86   24.09   30.55 ; ...
          9.21   10.26   13.31   16.49   19.71   22.95   26.18   32.68 ; ...
         11.36   12.42   15.47   18.67   21.90   25.15   28.40   34.92 ; ...
         13.64   14.69   17.76   20.96   24.21   27.47   30.74   37.29 ; ...
         16.04   17.09   20.16   23.39   26.65   29.93   33.21   39.80 ; ...
         18.57   19.62   22.70   25.94   29.22   32.52   35.82   42.44 ; ...
         21.24   22.29   25.38   28.64   31.94   35.25   38.58   45.24 ; ...
         24.05   25.10   28.21   31.49   34.80   38.14   41.48   48.18 ; ...
         27.02   28.07   31.20   34.49   37.83   41.18   44.55   51.30 ; ...
         30.16   31.21   34.34   37.66   41.02   44.40   47.79   54.58 ; ...
         33.47   34.52   37.67   41.00   44.39   47.79   51.20   58.04 ; ...
         36.96   38.01   41.17   44.53   47.94   51.37   54.81   61.70 ; ...
         40.65   41.69   44.88   48.26   51.69   55.15   58.61   65.56 ; ...
         44.54   45.58   48.79   52.19   55.65   59.13   62.63   69.63 ; ...
         48.65   49.69   52.91   56.34   59.83   63.34   66.87   73.93 ; ...
         52.99   54.03   57.27   60.73   64.25   67.79   71.34   78.47 ; ...
         57.57   58.61   61.87   65.36   68.91   72.48   76.07   83.26 ; ...
         62.41   63.46   66.74   70.25   73.83   77.44   81.06   88.32 ; ...
         67.53   68.58   71.88   75.42   79.04   82.68   86.33   93.66 ; ...
         72.94   73.98   77.31   80.89   84.54   88.21   91.90   99.31 ; ...
         78.66   79.70   83.05   86.66   90.35   94.06   97.79  105.27 ; ...
         84.71   85.75   89.12   92.77   96.49  100.24  104.01  111.57 ; ...
         91.11   92.15   95.55   99.23  102.99  106.78  110.59  118.24 ; ...
         97.88   98.93  102.34  106.06  109.86  113.70  117.55  125.28 ; ...
        105.05  106.09  109.54  113.29  117.14  121.02  124.91  132.73 ; ...
        112.65  113.69  117.15  120.95  124.84  128.76  132.70  140.62 ; ...
        120.69  121.73  125.22  129.06  132.99  136.96  140.95  148.96 ; ...
        129.22  130.25  133.77  137.65  141.63  145.65  149.69  157.80 ; ...
        138.25  139.29  142.83  146.76  150.78  154.85  158.94  167.16 ; ...
        147.84  148.87  152.44  156.41  160.49  164.61  168.75  177.08 ; ...
        158.01  159.04  162.64  166.65  170.78  174.96  179.16  187.59 ; ...
        168.81  169.83  173.46  177.52  181.70  185.93  190.19  198.74 ; ...
        ];

    pColdEos = [V0 K0 KP0];
    pHotEos  = [Tdeb0 gam0 a b 1.0];
    
    coldEosFun = @VinetEos;
    debyeDerivsFun = @debyeTange;
    hotEosFun  = @(V,T,V0,T0,pHotEos)...
        (MieGrunDebyeHotEos(V,T,V0,T0,Natom,pHotEos,debyeDerivsFun));
    elecThermPressFun = [];

    [P,KT,Cv,gam] = calcPressThermAddEos(V(:),T(:),T0,pColdEos,pHotEos,...
        coldEosFun,hotEosFun,elecThermPressFun);
    PmodTbl = reshape(P,length(Vrow),[]);

    calcErr = PmodTbl-PTbl;
    realInd = find(~isnan(calcErr));

    verifyTrue(testCase,all(reshape(abs(calcErr(realInd)),[],1) < TOL));
end

