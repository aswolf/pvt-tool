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
    peos = [V0,100,4.5];
    verifyEqual(testCase,VinetEos(V0,peos),0);
end
function testVinetEos_K0(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 100;
    peos = [V0,K0,4.5];
    dVstep = .001*[-1,1];
    V = V0*(1+dVstep);
    P = VinetEos(V,peos);
    K0num = -V0*diff(P)/diff(V);

    calcErr = K0num/K0-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testVinetEos_KP0(testCase)
    TOL = 1e-5;

    V0 = 10;
    K0 = 100;
    KP0 = 4.5;
    peos = [V0,K0,KP0];
    dVstep = 1e-4*[-1:1:1];
    V = V0*(1+dVstep);
    P = VinetEos(V,peos);

    % Perform 2nd order Taylor expansion to determine KP0
    K0num = -V0*diff(P)/diff(V);
    designM = [.5*P(:).^2 P(:) ones(size(P(:)))];
    polyv = designM\V(:);
    KP0num = V0*polyv(1)/polyv(2)^2 -1;

    calcErr = KP0num/KP0-1;
    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
function testVinetEos_Keqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 200;
    peos = [V0,K0,4.5];

    V = V0*.7*(1+1e-4*[-2:2]);
    [P dE K] = VinetEos(V,peos);

    Knum = -V.*central_diff(P,V);
    calcErr = Knum(3)./K(3)-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end

function testVinetEos_KPeqn(testCase)
    TOL = 1e-5;
    V0 = 10;
    K0 = 200;
    peos = [V0,K0,4.5];

    V = V0*.7*(1+1e-4*[-2:2]);
    [P dE K KP] = VinetEos(V,peos);

    KPnum = central_diff(K,P);
    calcErr = KPnum(3)./KP(3)-1;

    verifyTrue(testCase,all(abs(calcErr) < TOL));
end
