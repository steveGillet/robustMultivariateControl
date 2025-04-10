s = tf('s');
G11 = 2/(s+1);
G12 = 1/((s+1)*(s+2));
G21 = 1/((s+1)*(s+2));
G22 = 2/(s+2);

G = [G11 G12; G21 G22];

K11 = 2/(s+2);
K22 = 1/(s+3);
K = [K11 0; 0 K22];

L = G*K;
I = eye(2);
So = inv(I + L);
To = L*So;

figure;
[svG, w] = sigma(G);
semilogx(w, svG(2,:), w, svG(1,:));
title('Singular Values of G');
legend('Minimum Singular Values', 'Maximum Singular Values');
figure;
[svTo, w] = sigma(To);
semilogx(w, svTo(2,:), w, svTo(1,:));
title('Singular Values of To');
legend('Minimum Singular Values', 'Maximum Singular Values');
figure;
[svSo, w] = sigma(So);
semilogx(w, svSo(2,:), w, svSo(1,:));
title('Singular Values of So');
legend('Minimum Singular Values', 'Maximum Singular Values');

Mout = -To;
Min = -K*So*G;
gammaOut = norm(Mout, inf);
gammaIn = norm(Min, inf);

[svMout, wMout] = sigma(Mout);
[svMoutMax, iOut] = max(svMout(1,:));
freqOut = wMout(iOut);

disp(1/gammaOut);
disp(freqOut);

[svMin, wMin] = sigma(Min);
[svMinMax, iIn] = max(svMin(1,:));
freqIn = wMin(iIn);

disp(1/gammaIn);
disp(freqIn);

tau = 50;
Gscalar = 1/(tau*s + 1);
Gmatrix = [-87.8 1.4; -108.2 -1.4];
G = Gscalar*Gmatrix;

Kscalar = (tau * s + 1)/s;
Kmatrix = [-0.0015 0; 0 -0.075];
K = Kscalar*Kmatrix;

W1 = makeweight(0.1, 10, 1.2); 
W2 = makeweight(0.05, 50, 1.5);

figure;
opts = bodeoptions;
opts.MagUnits = 'abs'; % Absolute units instead of dB
opts.FreqUnits = 'rad/s';
bodeplot(W1, 'b-', W2, 'r--', opts, {1e-2, 1e3}); % Frequency range: 0.01 to 1000 rad/s
grid on;
title('Magnitude of Weighting Functions');
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
legend('W_1(s)', 'W_2(s)', 'Location', 'best');

W = blkdiag(W1, W2);
wMout = W * Mout;
wMin = W * Min;

normWmOut = norm(wMout, inf);
normWmIn = norm(wMin, inf);

betaOut = 1/normWmOut;
betaIn = 1/normWmIn;

disp(normWmOut);
disp(normWmIn);
disp(betaOut);
disp(betaIn);

figure;
subplot(2,1,1);
p = sigmaoptions('cstprefs');
p.MagUnits = 'abs';
sigmaplot(wMout, p);
hold on;
title('Singular Values of W(s)M(s) - Multiplicative Output Uncertainty');
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
sigmaplot(wMin, p);
hold on;
title('Singular Values of W(s)M(s) - Inverse Multiplicative Input Uncertainty');
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
grid on;
