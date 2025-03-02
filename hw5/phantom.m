A = [ -0.8 -0.0006 -12 0; 0 -0.014 -16.64 -32.2; 1 -0.0001 -1.5 0; 1 0 0 0];
B = [ -19 -2.5; -0.66 -0.5; -0.16 -0.6; 0 0];
C = [1 0 0 0; 0 0 1 0];
D = [0 0; 0 0];

phantomSS = ss(A,B,C,D);
[sv,wout] = sigma(phantomSS, {0,10});

figure;
plot(wout, sv, 'Color', [0.545, 0, 0]);
grid on;
title('Singular Value Bode Plot');
xlabel('Frequency [rad/s]');
ylabel('Absolute Gain [absolute units]');
hold on;

[num, den] = ss2tf(A,B,C,D,1);
tf11 = tf(num(1,:), den);
tf12 = tf(num(2,:), den);
[num, den] = ss2tf(A,B,C,D,2);
tf21 = tf(num(1,:), den);
tf22 = tf(num(2,:), den);

[mag11, phase11, wout11] = bode(tf11, {0,10});
[mag12, phase12, wout12] = bode(tf12, {0,10});
[mag21, phase21, wout21] = bode(tf21, {0,10});
[mag22, phase22, wout22] = bode(tf22, {0,10});
plot(wout11, squeeze(mag11), 'Color', [1,0.65,0]);
plot(wout12, squeeze(mag12), 'Color', [1,0.65,0]);
plot(wout21, squeeze(mag21), 'Color', [1,0.65,0]);
plot(wout22, squeeze(mag22), 'Color', [1,0.65,0]);
legend('System Singular Values', '', 'Individual Transfer Function Magnitudes');