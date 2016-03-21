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

% beta = @(MolCamPerMolPP2B, totalCam) (MolCamPerMolPP2B*totalPP2B -totalCam)./(totalCam-totalPP2B);
% KD5 = @(beta, totalCam) totalCam./beta - totalPP2B./(1+beta);
% kd5 = KD5(beta(y, t), t);
% 
% KD5 = @(MolCamPerMolPP2B, beta) totalPP2B*MolCamPerMolPP2B./(beta.*(1+beta));
% kd5_2 = KD5(y, beta(y, t));


KD5 = @(totalCam, MolCamPerMolPP2B) (totalCam - MolCamPerMolPP2B.*totalPP2B).*(1-MolCamPerMolPP2B)./MolCamPerMolPP2B;
kd5 = KD5(t, y);
kd5 = sort(kd5);
startmu = median(kd5);


MolCaMPerMolPP2B = @(KD5, totalCaM) (-(totalCaM.*(KD5 - totalPP2B + totalCaM - (KD5.^2 ...
    + 2*KD5.*totalPP2B + 2*KD5.*totalCaM + totalPP2B^2 ...
    - 2*totalPP2B.*totalCaM + totalCaM.^2).^(1/2)))./(totalPP2B.*(KD5 ...
    + totalPP2B - totalCaM + (KD5.^2 + 2*KD5.*totalPP2B + 2*KD5.*totalCaM ...
    + totalPP2B^2 - 2*totalPP2B*totalCaM + totalCaM.^2).^(1/2))));

f = @(x, mu, y) MolCaMPerMolPP2B(mu,x) - y;
f_res = @(mu) f(t, mu, y);
lb = kd5(1);
ub = kd5(end);
[x,resnorm,residual]  = lsqnonlin(f_res, startmu, lb, ub)
startmu = 26370
[x,resnorm,residual]  = lsqnonlin(f_res, startmu)


%Genom att ändra i1 och i2 styr man vilka som plottas
i1 = 5;
i2 = 8;

x1 = t(i1);
y1 = y(i1);


x2 = t(i2);
y2 = y(i2);



muguess = 26370;
mu = linspace(muguess*0.8, muguess*1.2);%26370
mu = linspace(muguess*0, muguess*2)
plot(f(x1, mu, y1), f(x2, mu, y2))
str = ['f_', num2str(i2),' as a funktion of f_', num2str(i1)];
title(str)
disp('press space to continue')
pause;

imax = length(t);
j = length(mu);
fval = zeros(imax, j);


%============================================%
%Vi testade att göra det för alla x-värden:
for i = 1:imax
    %fval(i, :) = f(t(i),mu,y(i));
    %ovanstående rad kan användas för att se att samtliga funktioner är
    %monotona.
    
    fval(i, :) = abs(f(t(i),mu,y(i)));
end

[val, muIndex] = min(max(abs(fval)));
%Detta tar inte hänsyn till om någon punkt är mycket dålig, alltså inte
%borde få så mycket utrymme. 
optMu = mu(muIndex)

axismax = 0.1;
hold off;  
for i = 1:imax
    plot(mu, fval(i, :))
    hold on;
    plot(mu(muIndex), fval(i, muIndex), 'o');
    plot(mu(muIndex)*ones(20,1), linspace(0, axismax, 20), ':')
    %hold off;
    axis([muguess*0 muguess*2 0 axismax])
    str = ['i = ', num2str(i)];
    title(str)
    pause
end

muIndex_left = muIndex - 1;
muIndex_right = muIndex + 1;

[~, i_max_left] = max(fval(:, muIndex_left));
[~, i_max_right] = max(fval(:, muIndex_right));

plot(mu, fval(i_max_left, :), 'r');
plot(mu, fval(i_max_right, :), 'r');

disp('Det ter sig som att de värsta värdena är xn, där n är:')
n1 = i_max_left
n2 = i_max_right

disp('Det optimala värdet på mu styrs främst av dessa, förutsatt att vi vill minimera det största felet')
disp('Det optimala värdet på mu blir då')
optMu


