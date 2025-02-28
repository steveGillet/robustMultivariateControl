% Problem 1

A = [-5  1  0;    % 3x3 matrix A
      0  1  1;
      1  1  1];

B2 = [0  0;       % 3x2 matrix B2
       0  1;
       1  0];

B1 = [0.5;         % 3x1 matrix B1
       0;
       0.3];

C1 = [1  0  0;    % 2x3 matrix C1
       0  2  1];

C2 = [1 0 0;
       0 0 1];

D21 = [0; 0];

nStates = size(A, 1);
nInputs = size(B2, 2);
nOutputs = size(C1, 1);

setlmis([])

P = lmivar(2, [nStates nStates]);
W = lmivar(2, [nStates nInputs]);
gamma = 0.1;

lmiterm([-1 1 1 P],1,1);

lmiterm([2 1 1 P],1,A,'s');
lmiterm([2 1 1 W],1, C2, 's');
lmiterm([2 1 2 P],1,B1);
lmiterm([2 1 3 0],C1');
lmiterm([2 2 2 0],-gamma);
lmiterm([2 2 3 0],0);
lmiterm([2 3 3 0],-gamma);

lmiH2obs = getlmis;

[tmin, xfeas] = feasp(lmiH2obs);
P = dec2mat(lmiH2obs,xfeas,P);
W = dec2mat(lmiH2obs,xfeas,W);

L = inv(P)*W;

hinfObsSys = ss(A+L*C2, B1, C1, 0);

[z, tOut, x] = impulse(hinfObsSys);

figure;

subplot(1,3,1);
plot(tOut, x(:,1), 'LineWidth', 1.5, 'Color', [255/255 127/255 102/255]);
grid on;
title('e1 Response to Impulse w');
xlabel('Time [s]');
ylabel('e1 Magnitude');

subplot(1,3,2);
plot(tOut, x(:,2), 'r-');
grid on;
title('e2 Response to Impulse w');
xlabel('Time [s]');
ylabel('e2 Magnitude');

subplot(1,3,3);
plot(tOut, x(:,3), 'y-');
grid on;
title('e3 Response to Impulse w');
xlabel('Time [s]');
ylabel('e3 Magnitude');
