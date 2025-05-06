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

nOut = size(C,1);
nIn = size(B,2);

M = 2;
Aweight = 0.005;
omegaB = 1;

s = tf('s');
Wp_theta = (s/100 + 1)/(s + 1*0.01);  % Low-frequency gain ~100, bandwidth ~1 rad/s
Wp_psi = (s/100 + 1)/(s + 1*0.01);    % Adjust parameters later
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
