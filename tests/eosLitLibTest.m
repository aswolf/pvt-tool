function tests = eosLitLibTest
    tests = functiontests(localfunctions);
end
function testGetEos_MgOTange2009(testCase)
    eosMgOMod = getEos_MgOTange2009();
    V0 = eosMgOMod.pEos(1);

    Tcol = [300    500   1000   1500   2000   2500   3000   4000];

    Vrow = V0*[1.15 1.14 1.13 1.12 1.11 1.10 1.09 1.08 1.07 1.06 1.05 ...
        1.04 1.03 1.02 1.01 1.00 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 ...
        0.91 0.90 0.89 0.88 0.87 0.86 0.85 0.84 0.83 0.82 0.81 0.80 0.79 ...
        0.78 0.77 0.76 0.75 0.74 0.73 0.72 0.71 0.70 0.69 0.68 0.67 0.66 ...
        0.65];

    [TM,VM] = meshgrid(Tcol,Vrow);
    V = VM(:);
    T = TM(:);

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


    [P] = evalPressEos([],eosMgOMod,V,T);
    PmodTbl = reshape(P,length(Vrow),[]);

    calcErr = PmodTbl-PTbl;
    realInd = find(~isnan(calcErr));

    verifyTrue(testCase,all(reshape(abs(calcErr(realInd)),[],1) < TOL),...
        'Calculated Pressures must agree with Table vals within TOL');
end
function testGetEos_MgPvTange2012(testCase)
    eosMgPvMod = getEos_MgPvTange2012();
    V0 = eosMgPvMod.pEos(1);

    Tcol = [300 1000 2000 3000 4000];
    Vrow = V0*[1.00 0.98 0.96 0.94 0.92 0.90 0.88 0.86 0.84 0.82 0.80 0.78 0.76];
    [TM,VM] = meshgrid(Tcol,Vrow);
    T = TM(:);
    V = VM(:);

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

    [P] = evalPressEos([],eosMgPvMod,V,T);
    PmodTbl = reshape(P,length(Vrow),[]);

    calcErr = PmodTbl-PTbl;

    realInd = find(~isnan(calcErr));
    verifyTrue(testCase,all(reshape(abs(calcErr(realInd)),[],1) < TOL),...
        'Calculated Pressures must agree with Table vals within TOL');
end
function testGetPVTdata_MgPvTange2012(testCase)
    PVTdata = getPVTdata_MgPvTange2012();
    verifyTrue(testCase,true,'Must be able to initialize PVTdata without error');
end
