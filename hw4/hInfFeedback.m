% % Problem 1

% A = [-5  1  0;    % 3x3 matrix A
%       0  1  1;
%       1  1  1];

% B2 = [0  0;       % 3x2 matrix B2
%        0  1;
%        1  0];

% B1 = [0.5;         % 3x1 matrix B1
%        0;
%        0.3];

% C1 = [1  0  0;    % 2x3 matrix C1
%        0  2  1];
  
% Problem 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [-0.0869 0 0.039 -1 ; % x1 = delta beta (rad)
   -4.424 -1.184 0 0.335; % x2 = delta p (rad/s)
                0 1 0 0 ; % x3 = delta phi (rad)
  2.148 -0.021 0 -0.228]; % x4 = delta r (rad/s)
B1 = [0 0 0 0.288]'; % w = r_g (yaw rate gust, rad/s)
B2 = [0.0223 0.547 0 -1.169]'; % u1 = delta_r
C1 = eye(4); % C = I
D11 = zeros(4,1);
D12 = [0 0 0 0]';

nStates = size(A, 1);
nInputs = size(B2, 2);
nOutputs = size(C1, 1);

openSS = ss(A,B1,C1,0);
disp('Open Loop Hinf Norm: ');
disp(hinfnorm(openSS));

setlmis([])

X = lmivar(1, [nStates 1]);
W = lmivar(2, [nInputs nStates]);
gamma = 2.5;

lmiterm([-1 1 1 X],1,1);
lmiterm([2 1 1 X],A,1,'s');
lmiterm([2 1 1 W],B2,1,'s');
lmiterm([2 1 2 0],B1);
lmiterm([2 1 3 X],1,C1');
lmiterm([2 2 2 0],-gamma);
lmiterm([2 2 3 0],0);
lmiterm([2 3 3 0],-gamma);

lmiDstab = getlmis;

[tmin, xfeas] = feasp(lmiDstab);
X = dec2mat(lmiDstab,xfeas,X);
W = dec2mat(lmiDstab,xfeas,W);

disp('X: ');
disp(X);
disp('W: ');
disp(W);

K = W/X;
disp('K: ');
disp(K);

hinfSys = ss(A+B2*K, B1, C1, 0);
disp("Hinf Norm");
disp(hinfnorm(hinfSys));

disp('Eigenvalues A+B2*K: ');
disp(eig(A+B2*K));

[z, tOut, x] = impulse(hinfSys);

figure;

subplot(1,2,1);
plot(tOut, z(:,1), 'LineWidth', 1.5, 'Color', [255/255 127/255 102/255]);
grid on;
title('z1 Response to Impulse w');
xlabel('Time [s]');
ylabel('z1 Magnitude');

subplot(1,2,2);
plot(tOut, z(:,2), 'r-');
grid on;
title('z2 Response to Impulse w');
xlabel('Time [s]');
ylabel('z2 Magnitude');

[sv,wout] = sigma(hinfSys);

figure;
% plot(wout, 10.^(sv/20));
plot(wout, sv);
grid on;
title('Max Singular Value at Different Frequencies');
xlabel('Frequency [rad/s]');
ylabel('Absolute Gain [absolute units]');

u = K*x';

% figure;
% plot(tOut, u(1,:), 'b-', tOut, u(2,:), 'g-');
% legend('u1', 'u2');
% grid on;
% title('Control Signal Response to Impulse Disturbance');
% xlabel('Time [s]');
% ylabel('Control Signal Magnitude');

disp("u1 L2 Norm: ")
disp(norm(u(1,:),2));

% disp("u2 L2 Norm: ")
% disp(norm(u(2,:),2));

figure

subplot(1,4,1);
plot(tOut, x(:,1), 'LineWidth', 1.5, 'Color', [255/255 127/255 102/255]);
grid on;
title('x1 Response to Impulse w');
xlabel('Time [s]');
ylabel('x1 Magnitude');

subplot(1,4,2);
plot(tOut, x(:,2), 'r-');
grid on;
title('x2 Response to Impulse w');
xlabel('Time [s]');
ylabel('x2 Magnitude');

subplot(1,4,3);
plot(tOut, x(:,3), 'r-');
grid on;
title('x3 Response to Impulse w');
xlabel('Time [s]');
ylabel('x3 Magnitude');

subplot(1,4,4);
plot(tOut, x(:,3), 'r-');
grid on;
title('x4 Response to Impulse w');
xlabel('Time [s]');
ylabel('x4 Magnitude');

disp("Average Power in States: ");
Phinf = lyap(A+B2*K, B1 * B1');
averagePower = trace(Phinf);
disp(averagePower);