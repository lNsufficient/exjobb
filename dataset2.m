clear;
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
% MolCaMPerMolPP2B = @(KD5, totalCaM) (-(totalCaM.*(KD5 - totalPP2B + totalCaM - (KD5.^2 ...
%     + 2*KD5.*totalPP2B + 2*KD5.*totalCaM + totalPP2B^2 ...
%     - 2*totalPP2B.*totalCaM + totalCaM.^2).^(1/2)))/(totalPP2B*(KD5 ...
%     + totalPP2B - totalCaM + (KD5.^2 + 2*KD5*totalPP2B + 2*KD5*totalCaM ...
%     + totalPP2B^2 - 2*totalPP2B*totalCaM + totalCaM.^2).^(1/2))));
% 
MolCaMPerMolPP2B = @(KD5, totalCaM) (-(totalCaM.*(KD5 - totalPP2B + totalCaM - (KD5.^2 ...
    + 2*KD5.*totalPP2B + 2*KD5.*totalCaM + totalPP2B^2 ...
    - 2*totalPP2B.*totalCaM + totalCaM.^2).^(1/2)))./(totalPP2B.*(KD5 ...
    + totalPP2B - totalCaM + (KD5.^2 + 2*KD5.*totalPP2B + 2*KD5.*totalCaM ...
    + totalPP2B^2 - 2*totalPP2B*totalCaM + totalCaM.^2).^(1/2))));

nbrRuns = 4;

res = zeros(nbrRuns, 20);
mu = zeros(nbrRuns, 1);

gaussIndex = 1;
medianIndex = 2;
matlabIndex = 3;
startIndex = 4;

mu(startIndex) = 26370;



f = @(x, mu, y) MolCaMPerMolPP2B(mu,x) - y;

KD5 = @(totalCam, MolCamPerMolPP2B) (totalCam - MolCamPerMolPP2B.*totalPP2B).*(1-MolCamPerMolPP2B)./MolCamPerMolPP2B;
kd5 = KD5(t, y);
kd5 = sort(kd5);
mu(medianIndex) = median(kd5);

res(startIndex, :) = (f(t,mu(startIndex), y));
res(medianIndex, :) = (f(t,mu(medianIndex), y));

mu(gaussIndex) = gaussnewtonmod(MolCaMPerMolPP2B,t,y,mu(medianIndex),1e-4,1,1,1);
% f = @(alpha) sum((MolCaMPerMolPP2B(alpha,t) - y).^2);
% sign=f(mu(1))-f(26370);
% if sign < 0
%     disp('improvment')
% end


res(gaussIndex, :) = (f(t,mu(gaussIndex), y));



f_res = @(mu) f(t, mu, y);
lb = kd5(1);
ub = kd5(end);
[mu(matlabIndex), ~, res(matlabIndex,:)]  = lsqnonlin(f_res, mu(medianIndex), lb, ub);

%startmu = 26370
%[mu_matlab,resnorm,residual]  = lsqnonlin(f_res, startmu)
%plot(residual)


%26370
%x = linspace(min(t), max(t));

clf
plot(t, y, '.');
hold on;
%t = x;
plot(t, MolCaMPerMolPP2B(mu(matlabIndex),t));

plot(t, MolCaMPerMolPP2B(mu(gaussIndex),t), 'r');
plot(t, MolCaMPerMolPP2B(mu(medianIndex),t), 'm');
plot(t, MolCaMPerMolPP2B(mu(startIndex),t), ':k');

plot(x, MolCaMPerMolPP2B(mu(startIndex),x), ':k');
clf
semilogy(t, abs(res))

norms = sqrt(sum(abs(res).^2,2));

	