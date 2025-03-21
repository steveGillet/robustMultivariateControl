Lp = -1;
Lda = 30;
tau = 0.1;

A = [Lp 0 Lda; 1 0 0; 0 0 -1/tau];
B = [0; 0; 1/tau];
C = [0 1 0];
D = 0;
G = ss(A,B,C,D);

M = 2;
A = 0.005;
omegaB = .5;

W1 = tf([1/M omegaB], [1 omegaB*A]);
W2 = tf([100 10], [1 100]);

systemnames = 'G W1 W2';
inputvar = '[w; u]';
outputvar = '[W1; W2; w-G]';
input_to_G = '[u]';
input_to_W1 = '[w-G]';
input_to_W2 = '[u]';
sysoutname = 'P';
sysic;
P = minreal(ss(P));

m = size(B,2);
p = size(C,1);
[K, CL, gamma, info] = hinfsyn(P,p,m,'method','ric','Tolgam',1e-3,'DISPLAY','on');

Si = minreal(inv(eye(1)+K*G));
So = minreal(inv(eye(1)+G*K));
GSi = minreal(G*inv(eye(1)+K*G));

figure(1); clf;
sigma(G, 'r', K, 'y', G*K, 'm', {1e-2, 1e2});
grid on; legend('\sigma(G)','\sigma(K)','\sigma(G*K)'); title('Open Loop Singular Values');

figure(2); clf;
sigma(So, 'r', inv(W1), 'y', {1e-3, 1e2});
grid on; legend('\sigma(So)','\sigma(W1)'); title('Closed Loop Singular Values');

figure(3); clf;
sigma(K*So, 'r', inv(W2), 'y', {1e-3, 1e2});
grid on; legend('\sigma(KSo)','\sigma(W2)'); title('Closed Loop Singular Values');

figure(4); clf;
prefilter = dcgain(K);
step(prefilter*CL); grid on; axis([0 3 0 1.2]);
hinfnorm(prefilter*CL);
title('Closed Loop Step Response');
