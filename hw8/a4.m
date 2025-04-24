Ybeta = -110.94;
u0 = 1126;
Yr = 0;
Nbeta = 15.16;
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
A = 0.005;
omegaB = 1;

Wp = tf([1/M omegaB], [1 omegaB*A]);

systemnames = 'G Wp';
inputvar = '[deltaa; deltar]';
outputvar = '[beta; r]';
input_to_G = '[deltaa; deltar]';
input_to_Wp = '[deltaa; deltar]';
sysoutname = 'P';
sysic;
P = minreal(ss(P));

m = size(B,2);
p = size(C,1);
[K, CL, gamma, info] = hinfsyn(P,p,m,'method','ric','Tolgam',1e-3,'DISPLAY','on');
