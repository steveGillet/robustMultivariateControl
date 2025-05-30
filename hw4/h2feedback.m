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

setlmis([])

X = lmivar(1, [nStates 1]);
W = lmivar(2, [nInputs nStates]);
Z = lmivar(1, [nStates 1]);
gamma = 1.14;

lmiterm([-1 1 1 X],1,1);

lmiterm([2 1 1 X],A,1,'s');
lmiterm([2 1 1 W],B2,1,'s');
lmiterm([2 1 1 0],B1*B1');

lmiterm([3 1 1 Z],-1,1); 
lmiterm([3 1 2 X],C1,1);
lmiterm([3 1 2 W],D12,1);
lmiterm([3 2 2 X],-1,1);

lmiterm([4 1 1 Z],[1 0 0 0], [1 0 0 0]');
lmiterm([4 1 1 Z],[0 1 0 0], [0 1 0 0]');
lmiterm([4 1 1 Z],[0 0 1 0], [0 0 1 0]');
lmiterm([4 1 1 Z],[0 0 0 1], [0 0 0 1]');
lmiterm([-4 1 1 0],gamma);

h2lmi = getlmis;

[tmin, xfeas] = feasp(h2lmi);
X = dec2mat(h2lmi,xfeas,X);
W = dec2mat(h2lmi,xfeas,W);
Z = dec2mat(h2lmi,xfeas,Z);

K = W/X;

h2sys = ss(A+B2*K,B1,C1,D11);

[z, tOut, x] = impulse(h2sys);

u = K*x';

disp("u1 L2 Norm: ")
disp(norm(u(1,:),2));

% evalsys = evallmi(h2lmi,xfeas);
% [lhs,rhs] = showlmi(evalsys,4);

disp('H2 norm');
disp(norm(h2sys, 2));

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

[sv,wout] = sigma(h2sys);

figure;
% plot(wout, 10.^(sv/20));
plot(wout, sv);
grid on;
title('Max Singular Value at Different Frequencies');
xlabel('Frequency [rad/s]');
ylabel('Absolute Gain [absolute units]');

disp("Average Power in States: ");
Ph2 = lyap(A+B2*K, B1 * B1');
averagePower = trace(Ph2);
disp(averagePower);