% Testing Gauss-Newton
%[Ca]	mol Ca/mol CaM
Matrix=[457.088	0.166;
831.763	0.342;
1122.018	0.528;
1513.561	0.753;
1949.844	0.909;
2398.832	1.115;
2818.382	1.311;
3467.368	1.477;
4168.693	1.682;
5128.613	1.839;
5754.399	2.015;
7585.775	2.142;
8317.637	2.328;
11481.536	2.426;
14791.083	2.69;
15135.612	2.857;
26302.679	3.111;
26915.348	3.287;
39810.717	3.395;
40738.027	3.551;
53703.179	3.60;
70794.578	3.649;
93325.430	3.855;
151356.124	3.796;
194984.459	4.021;
257039.578	3.894;
346736.850	3.953;
524807.460	4.001;
676082.975	4.06;]
t = (Matrix(:,1));
y = (Matrix(:,2));
MolCaPerMolCaM = @(KD, Ca)(Ca.*(4*Ca.^3 + 3*KD(4)*Ca.^2 + 2*KD(3)*KD(4)*Ca...
    + KD(2)*KD(3)*KD(4)))./(Ca.^4 + KD(4)*Ca.^3 + KD(3)*KD(4)*Ca.^2 +...
    KD(2)*KD(3)*KD(4)*Ca + KD(1)*KD(2)*KD(3)*KD(4));
    x = gaussnewton(MolCaPerMolCaM,t,y,[4646, 1124, 22959, 15836]',1e-6,0,1,1)
    f = @(alpha) sum((MolCaPerMolCaM(alpha,t) - y).^2);
    f(x)-f([4646, 1124, 22959, 15836]')
  