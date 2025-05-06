% Parameters for the 2-wheeled self-balancing robot with differential drive
M = 1.5;        % Body mass (kg)
m = 0.3;        % Total wheel mass (kg)
l = 0.15;       % Distance from wheel axle to center of mass (m)
R = 0.075;      % Wheel radius (m)
d = 0.15;       % Wheel separation (m)
I_psi = 0.0141; % Moment of inertia about vertical axis (kg·m^2)
g = 9.81;       % Gravitational acceleration (m/s^2)

% State-space matrices
A = [0, 1, 0, 0, 0;
     -(M + m) * g / (M * l), 0, 0, 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 0;
     m * g / M, 0, 0, 0, 0];

B = [0, 0;
     -2 / (M * l * R), 0;
     0, 0;
     0, d / (I_psi * R);
     2 / (M * R), 0];

% Define the state-space system
% Output matrix C assumes all states are measurable, D = 0
C = [1 0 0 0 0; 0 0 1 0 0];
D = zeros(2, 2);
G = ss(A, B, C, D);

% Check controllability and observability
disp('Controllability rank:');
disp(rank(ctrb(G.A, G.B)));
disp('Observability rank:');
disp(rank(obsv(G.A, G.C)));


nOut = size(C,1);
nIn = size(B,2);

M = 100;
Aweight = 0.01;
omegaB = 1;

s = tf('s');
Wp_theta = (s/M + omegaB)/(s + omegaB*Aweight);
Wp_psi = (s/M + omegaB)/(s + omegaB*Aweight);  
Wp = blkdiag(Wp_theta, Wp_psi);
Wp = ss(Wp);

% Display the system
disp('State-space system:');
disp(G);

% Optional: Simulate open-loop response (unstable without control)
t = 0:0.01:2;  % Time vector (2 seconds)
u = zeros(length(t), 2);  % No input (open-loop)
x0 = [0.1; 0; 0; 0; 0];  % Initial condition: 0.1 rad tilt
[y, t, x] = lsim(G, u, t, x0);
clf;figure(1);
subplot(3,1,1);
plot(t, x(:,1));
title('Tilt Angle \theta (rad)');
subplot(3,1,2);
plot(t, x(:,3));
title('Yaw Angle \psi (rad)');
subplot(3,1,3);
plot(t, x(:,5));
title('Forward Velocity v (m/s)');


systemnames = 'G Wp';
inputvar = '[w(2); u(2)]';
outputvar = '[Wp; -G-w]';
input_to_G = '[u]';
input_to_Wp = '[w + G]';
sysoutname = 'P';
sysic;
P = minreal(ss(P));


[K, CL, gamma, info] = hinfsyn(P,nOut,nIn,'method','ric','Tolgam',1e-3,'DISPLAY','on');

figure(2);
sigma(G);
title('Open-loop Transfer Function G Sigma Plot');

figure(3);
sigma(G*K);
title('Open-loop Transfer Function G*K Sigma Plot');

[Kmix, gammaMix] = mixsyn(G, Wp, [], []);

So = inv(eye(nOut)+G*K);
To = minreal(G*K*inv(eye(nOut)+G*K));

figure(4);

sigma(So);
title('Sensitivity Function So Sigma Plot');
figure(5);
sigma(To);
title('Complementary Sensitivity Function To Sigma Plot');
figure(6);
sigma(Wp);
title('Weighting Function Wp Sigma Plot');

CL = feedback(G*K, eye(nOut));

% Plot step response for theta_r
figure(7);
step(CL(1,1), 'b');  % theta response to theta_r
hold on;
step(CL(2,1), 'r');  % psi response to theta_r (cross-coupling effect)
title('Response to Step in \theta_r');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\theta', '\psi');
hold off;

% Plot step response for psi_r
figure(8);
step(CL(2,2), 'b');  % psi response to psi_r
hold on;
step(CL(1,2), 'r');  % theta response to psi_r (cross-coupling effect)
title('Response to Step in \psi_r');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\psi', '\theta');
hold off;

S = feedback(eye(nOut), G*K);

% Plot disturbance response for d_theta
figure(9);
step(S(1,1), 'b');  % theta response to d_theta
hold on;
step(S(2,1), 'r');  % psi response to d_theta
title('Response to Disturbance d_\theta');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\theta', '\psi');
hold off;

% Plot disturbance response for d_psi
figure(10);
step(S(2,2), 'b');  % psi response to d_psi
hold on;
step(S(1,2), 'r');  % theta response to d_psi
title('Response to Disturbance d_\psi');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\psi', '\theta');
hold off;

Wi = (s+0.2)/(0.5*s+1)*eye(nIn);

systemnames = 'G Wp Wi';
inputvar = '[ydel(2); w(2); u(2)]';
outputvar = '[Wi; Wp; -G-w]';
input_to_G = '[u + ydel]';
input_to_Wp = '[G+w]';
input_to_Wi = '[u]';
sysoutname = 'P';
sysic;
P = minreal(ss(P));


d0 = 1;
D = append(d0, d0, tf(eye(2)), tf(eye(2)));
[K, CL, gamma(1), info] = hinfsyn(D*P*inv(D),nOut,nIn,'method','ric','Tolgam',1e-3,'DISPLAY','on');

omega = logspace(-3,3);
blk = [1 0; 1 0];

N = frd(lft(P,K),omega);
[muRS, infoRS] = mussv(N(1:2,1:2), blk);

bodeOptions = bodeoptions;
bodeOptions.PhaseVisible = 'off';
bodeOptions.XLim = [1e-3 1e3];
bodeOptions.MagUnits = 'abs';

figure(11); clf;
bodeplot(muRS, bodeOptions);
hold on; grid on;
title('Robust Stability Plot');

blk = [1 0; 1 0; 2 2];
[muRP, infoRP] = mussv(N, blk);
figure(12);
bodeplot(muRP, bodeOptions);
hold on; grid on;
title('Robust Performance Plot');

muRSpeak = max(sigma(muRS));
disp('Robust Stability Margin:');
disp(1/muRSpeak);

muRPpeak = max(sigma(muRP));
disp('Robust Performance Margin:');
disp(1/muRPpeak);

CL = feedback(G*K, eye(nOut));

% Plot step response for theta_r
figure(13);
step(CL(1,1), 'b');  % theta response to theta_r
hold on;
step(CL(2,1), 'r');  % psi response to theta_r (cross-coupling effect)
title('Response to Step in \theta_r');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\theta', '\psi');
hold off;

% Plot step response for psi_r
figure(14);
step(CL(2,2), 'b');  % psi response to psi_r
hold on;
step(CL(1,2), 'r');  % theta response to psi_r (cross-coupling effect)
title('Response to Step in \psi_r');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\psi', '\theta');
hold off;

S = feedback(eye(nOut), G*K);

% Plot disturbance response for d_theta
figure(15);
step(S(1,1), 'b');  % theta response to d_theta
hold on;
step(S(2,1), 'r');  % psi response to d_theta
title('Response to Disturbance d_\theta');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\theta', '\psi');
hold off;

% Plot disturbance response for d_psi
figure(16);
step(S(2,2), 'b');  % psi response to d_psi
hold on;
step(S(1,2), 'r');  % theta response to d_psi
title('Response to Disturbance d_\psi');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\psi', '\theta');
hold off;


Delta = [ultidyn('D_1', [1 1]) 0 ; 0 ultidyn('D_2', [1 1])];
Punc = lft(Delta,P);
[Kauto,clp,dkinfo] = musyn(Punc, nOut, nIn);

N = frd(lft(P,Kauto), omega);
[muAuto, infoAuto] = mussv(N, blk);
muRP = norm(muAuto(1,1),inf,1e-6);

disp('Robust Performance Margin for D-K Iterated Controller:');
disp(1/muRP);

figure(17); clf;
bodeplot(muAuto(1,1), bodeOptions); grid on;
title('Robust Performance Plot for D-K Iterated Controller');

So = inv(eye(nOut)+G*Kauto);
To = G*Kauto*inv(eye(nOut)+G*Kauto);

figure(18);
sigma(G,G*Kauto);
legend('G','G*K');
title('Open Loop Singular Values');

figure(19);
sigma(So, To, inv(Wp));
legend('So','To','1/Wp');
title('Closed Loop Singular Values');

CL = feedback(G*Kauto, eye(nOut));

% Plot step response for theta_r
figure(20);
step(CL(1,1), 'b');  % theta response to theta_r
hold on;
step(CL(2,1), 'r');  % psi response to theta_r (cross-coupling effect)
title('Response to Step in \theta_r');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\theta', '\psi');
hold off;

% Plot step response for psi_r
figure(21);
step(CL(2,2), 'b');  % psi response to psi_r
hold on;
step(CL(1,2), 'r');  % theta response to psi_r (cross-coupling effect)
title('Response to Step in \psi_r');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\psi', '\theta');
hold off;

S = feedback(eye(nOut), G*Kauto);

% Plot disturbance response for d_theta
figure(22);
step(S(1,1), 'b');  % theta response to d_theta
hold on;
step(S(2,1), 'r');  % psi response to d_theta
title('Response to Disturbance d_\theta');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\theta', '\psi');
hold off;

% Plot disturbance response for d_psi
figure(23);
step(S(2,2), 'b');  % psi response to d_psi
hold on;
step(S(1,2), 'r');  % theta response to d_psi
title('Response to Disturbance d_\psi');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\psi', '\theta');
hold off;

figure(24); clf;
subplot(2,2,1), step(G(1,1), To(1,1)); grid on;
subplot(2,2,2), step(G(2,2), To(2,2)); grid on;
subplot(2,2,3), step(So(1,1));
subplot(2,2,4), step(So(2,2));
