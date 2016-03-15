% Testing Gauss-Newton
%log[CaM]	[MolCaMPerMolPP2B]
totalPP2B = 100;
Matrix=[0.0     0.0
1.43	0.0;
1.83	0.0;
2.20	0.0;
2.61	0.027;
2.91	0.038;
3.22	0.085;
3.42	0.120;
3.73	0.158;
3.90	0.225;
4.08	0.304;
4.21	0.365;
4.38	0.473;
4.55	0.543;
4.69	0.669;
4.80	0.724;
4.93	0.800;
5.10	0.823;
5.30	0.885;
5.36	0.902];
t = 10.^(Matrix(:,1))';
y = (Matrix(:,2))';
MolCaMPerMolPP2B = @(KD5, totalCaM) (-(totalCaM.*(KD5 - totalPP2B + totalCaM - (KD5.^2 ...
    + 2*KD5.*totalPP2B + 2*KD5.*totalCaM + totalPP2B^2 ...
    - 2*totalPP2B.*totalCaM + totalCaM.^2).^(1/2)))/(totalPP2B*(KD5 ...
    + totalPP2B - totalCaM + (KD5.^2 + 2*KD5*totalPP2B + 2*KD5*totalCaM ...
    + totalPP2B^2 - 2*totalPP2B*totalCaM + totalCaM.^2).^(1/2))));
x = gaussnewtonmod(MolCaMPerMolPP2B,t,y,26370,1e-4,1,1,1)
f = @(alpha) sum((MolCaMPerMolPP2B(alpha,t) - y).^2);
sign=f(x)-f(26370)
if sign < 0
    disp('improvment')
end

%26370