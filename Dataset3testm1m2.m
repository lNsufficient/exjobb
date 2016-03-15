% Testing Gauss-Newton
%log[Ca](nM)	ActivePP2BPercentage
Matrix1= [2.65	10.47;
2.83	36.99;
3.01	67.89;
3.17	88.41;
3.49	97.74];
%log[Ca](nM)	ActivePP2BPercentage
Matrix2= [2.65	5.55;
2.83	11.31;
3.01	29.37;
3.17	53.43;
3.34	86.24;
3.49	90.64;
3.75	96.96;
3.99	99.72];

t = 10.^(Matrix1(:,1));
y = (Matrix1(:,2));

totalPP2B = 3; % Ska denna vara samma som innan?
totalCaM = 300; % Vad ska denna vara?

ActivePP2BPercentage = @(m1_7, m2_7, Ca) (200*(Ca.^4./(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7) ...
    - (((2*Ca.^4)./(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7) ...
    + (Ca.^4.*(totalPP2B - totalCaM))./(Ca.^4 + Ca.^3 + Ca.^2 + m2_7.*Ca + m1_7.*m2_7)).^2 ...
    + (8*Ca.^8*totalCaM)./((Ca.^4 + Ca.^3 + Ca.^2 + m2_7.*Ca + m1_7.*m2_7) ...
    .*(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7))).^(1/2)/2 ...
    + (Ca.^4*(totalPP2B - totalCaM))./(2*(Ca.^4 + Ca.^3 + Ca.^2 + m2_7.*Ca + m1_7.*m2_7)))) ...
    ./((((2*Ca.^4)./(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7) ...
    - (((2*Ca.^4)./(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7) ...
    + (Ca.^4*(totalPP2B - totalCaM))./(Ca.^4 + Ca.^3 + Ca.^2 + m2_7.*Ca + m1_7.*m2_7)).^2 ...
    + (8*Ca.^8*totalCaM)./((Ca.^4 + Ca.^3 + Ca.^2 + m2_7.*Ca + m1_7.*m2_7) ...
    .*(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7))).^(1/2) ...
    + (Ca.^4*(totalPP2B - totalCaM))./(Ca.^4 + Ca.^3 + Ca.^2 + m2_7.*Ca + m1_7.*m2_7)) ...
    .*(2*Ca.^4 + 2*Ca.^3 + 2*Ca.^2 + 2*m2_7.*Ca + m1_7.*m2_7))./(2*Ca.^4) - 2);


%[KD*CaM_Ca3*Ca] = KD4
%[KD*CaM_Ca2*Ca] = KD3
%[KD*CaM_Ca1*Ca] = KD2
%[KD*CaM*Ca] = KD1
%[KD*CaM_Ca4*PP2B] = KD9
%[KD*CaM_Ca3*PP2B] = KD8
%[KD*CaM_Ca2*PP2B] = KD7
%[KD*CaM_Ca1*PP2B] = KD6
%[KD*CaM*PP2B] = KD5
%[KD*PP2B_CaM_Ca3*Ca] = KD13 
%[KD*PP2B_CaM_Ca2*Ca] = KD12
%[KD*PP2B_CaM_Ca1*Ca] = KD11
%[KD*PP2B_CaM*Ca] = KD10

KD=[4646, 1124, 22959, 15836, 26370, 3405.5, 181.79, 6.3344, 0.028, 600, 60, 800, 70]';
Ca = 1000;

m1_7opt = @(KD, Ca) (2*Ca.*(Ca.^4*KD(5)*KD(6)*KD(7)*KD(8) + Ca.^3*KD(4)*KD(5)*KD(6)*KD(7)*KD(9) ...
    - Ca.^4*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) + Ca.^2*KD(3)*KD(4)*KD(5)*KD(6)*KD(8)*KD(9) ... 
    - Ca.^3*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) - Ca.^2*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) ...
    + Ca*KD(2)*KD(3)*KD(4)*KD(5)*KD(7)*KD(8)*KD(9) + KD(1)*KD(2)*KD(3)*KD(4)*KD(6)*KD(7)*KD(8)*KD(9) ...
    - KD(1)*KD(2)*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) - Ca*KD(2)*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9))) ...
    ./(Ca.^2*KD(5)*KD(6)*KD(7)*KD(8) + Ca.^3*KD(5)*KD(6)*KD(7)*KD(8) - Ca.^4*KD(5)*KD(6)*KD(7)*KD(8) ...
    - 2*Ca.^3*KD(4)*KD(5)*KD(6)*KD(7)*KD(9) + Ca.^4*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) ...
    - 2*Ca.^2*KD(3)*KD(4)*KD(5)*KD(6)*KD(8)*KD(9) + Ca.^3*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) ...
    + Ca.^2*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) - 2*Ca*KD(2)*KD(3)*KD(4)*KD(5)*KD(7)*KD(8)*KD(9) ...
    - 2*KD(1)*KD(2)*KD(3)*KD(4)*KD(6)*KD(7)*KD(8)*KD(9) + KD(1)*KD(2)*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) ...
    + Ca*KD(2)*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9));

m2_7opt =  @(KD, Ca) -(Ca.^2*KD(5)*KD(6)*KD(7)*KD(8) + Ca.^3*KD(5)*KD(6)*KD(7)*KD(8) ...
    - Ca.^4*KD(5)*KD(6)*KD(7)*KD(8) - 2*Ca.^3*KD(4)*KD(5)*KD(6)*KD(7)*KD(9) ...
    + Ca.^4*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) - 2*Ca.^2*KD(3)*KD(4)*KD(5)*KD(6)*KD(8)*KD(9) ...
    + Ca.^3*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) + Ca.^2*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) ...
    - 2*Ca*KD(2)*KD(3)*KD(4)*KD(5)*KD(7)*KD(8)*KD(9) - 2*KD(1)*KD(2)*KD(3)*KD(4)*KD(6)*KD(7)*KD(8)*KD(9) ...
    + KD(1)*KD(2)*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9) ...
    + Ca*KD(2)*KD(3)*KD(4)*KD(5)*KD(6)*KD(7)*KD(8)*KD(9))./(Ca*KD(5)*KD(6)*KD(7)*KD(8));

%m1_7 = @(KD, Ca) arrayfun(m1_7opt, KD.*ones(size(Ca)), Ca);

%m1 = m1_7opt(KD, t)
%m2 = m2_7opt(KD, t)

testf = @(KD, Ca) (ActivePP2BPercentage(m1_7opt(KD, Ca), m2_7opt(KD, Ca), Ca));
%testf = @(M, Ca) arrayfun(ActivePP2BPercentage,M(1)*ones(size(Ca)), M(2)*ones(size(Ca)), Ca);
%testf = @(M, Ca) ActivePP2BPercentage(M(1), M(2), Ca.^1)
startguess = KD;
x = gaussnewton(testf,t,y,startguess,1e-4,0,0,0)
f = @(alpha) sum((testf(alpha,t) - y).^2);
sign=f(x)-f(startguess)
if sign < 0
    disp('improvment')
end

xplot = linspace(0, max(t), 100);
plot(t, y, 'x');
hold on;
plot(xplot, testf(x, xplot))
plot(xplot, testf(KD, xplot))