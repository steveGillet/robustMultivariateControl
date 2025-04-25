Ybeta = -110.94;
Ybeta = ureal('Ybeta', Ybeta, 'Percentage', 20);
u0 = 1126*0.4;
Yr = 0;
Nbeta = 15.16;
Nbeta = ureal('Nbeta', Nbeta, 'Percentage', 30);
Nr = -0.639;
Ydeltar = 19.65;
Ndeltaa = 0.334;
Ndeltar = -6.732;

A = [Ybeta/u0 -(1-Yr/u0); Nbeta Nr];
B = [0 Ydeltar/u0; Ndeltaa Ndeltar];
C = eye(size(B,1));
D = 0;

G = ss(A,B,C,D);

M = 2;
Aweight = 0.005;
omegaB = 1;

s = tf('s');
Wp = tf([s/M+omegaB], [s+omegaB*Aweight])*eye(size(C,1));

systemnames = 'G Wp';
inputvar = '[u(2); w(2)]';
outputvar = '[Wp; -G-w]';
input_to_G = '[u]';
input_to_Wp = '[G+w]';
sysoutname = 'P';
sysic;
Pnominal = minreal(ss(P));

m = size(B,2);
p = size(C,1);
[K, CL, gamma, info] = hinfsyn(Pnominal,p,m,'method','ric','Tolgam',1e-3,'DISPLAY','on');

So = inv(eye(2)+G*K);
figure(1); clf;
sigma(So, 'm', inv(Wp), 'r--', {1e-3, 1e2});
grid on; legend('\sigma(So)','\sigma(1/Wp)'); title('\sigma(So) vs. |1/Wp(jw)|');

CL = lft(P,K);
[N, Delta, Blkstruct] = lftdata(CL);
szDelta = size(Delta);
N11 = N(1:szDelta(2), 1:szDelta(1));

omega = logspace(-3,2);
N11frd = frd(N11, omega);
mu = mussv(N11frd, Blkstruct);

bodeOpt = bodeoptions;
bodeOpt.PhaseVisible = 'off';
bodeOpt.XLim = [1e-1 1e1];
bodeOpt.MagUnits = 'abs';
figure(2); clf;
bodeplot(mu(1,1), 'bo', mu(1,2), 'r-', bodeOpt);
grid on;
xlabel('Frequency (rad/s)');
ylabel('Mu(N11) Upper and Lower Bounds');
title('Robust Stability Plot');

Nfrd = frd(N, omega);
Blkstruct = [1 0; 2 2];
muN = mussv(Nfrd, Blkstruct);
figure(3); clf;
bodeplot(muN(1,1), 'bo', muN(1,2), 'r-', bodeOpt);
grid on;
xlabel('Frequency (rad/s)');
ylabel('Mu(N) Upper and Lower Bounds');
title('Robust Performance Plot');
