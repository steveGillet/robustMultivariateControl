num = [40 1];
den = conv([1 1 4], [1 6]);
plant = tf(num, den);

bode(plant);

k = 125;
alpha = 10000;
tau = 0.25;
numC1 = [k*tau k*1];
denC1 = [alpha*tau 1];
C1 = tf(numC1, denC1);

k = 125;
alpha = 0.1;
tau = 0.5;
numC2 = [k*tau k];
denC2 = [alpha*tau 1];
C2 = tf(numC2, denC2);

openLoop = plant*C1*C2;
[~, ~, ~, Wc] = margin(openLoop);

figure;
bode(openLoop);
grid on;

ax = gca;
allAxes = findall(gcf, 'Type', 'Axes');

for i = 1:length(allAxes)
    axes(allAxes(i));
    hold on;
    xline(Wc, 'r-', 'Gain Crossover', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Alpha', 0.7);
    hold off;
end

axes(allAxes(1));
hold on;
xline(12.56, 'r--', 'Bandwidth Limit (2 Hz)', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Alpha', 0.7);
xline(3.14, 'g--', 'Tracking Error (0.5 Hz)', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Alpha', 0.7);
yline(6, 'g--', 'Open Loop Gain (6dB)', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Alpha', 0.7);
hold off;

axes(allAxes(2));
hold on;
yline(-135, 'r--', 'Phase Margin (45 degrees)', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'Alpha', 0.7);
hold off;

closedLoop = feedback(openLoop,1);

figure;
step(closedLoop);
grid on;
title('Closed Loop Step Response');
xlabel('Time [s]');
ylabel('Output \phi(t) [rad]');

t = 0:0.1:10;
omega = 2*pi*0.5;
r = sin(omega*t);
y = lsim(closedLoop, r, t);
e = r'-y;

figure;
plot(t,e, 'Color', [255/255 140/255 105/255]);
grid on;
title('Error Signal e(t) for r sin(2\pi\omega t) (\omega =0.5Hz)');
xlabel('Time [s]');
ylabel('Output \phi(t) [rad]');

figure;
plot(t,r, 'Color', [255/255 140/255 105/255]);
hold on;
plot(t,y, 'Color', [204/255 51/255 51/255]);
grid on;
title('r and y for r sin(2\pi\omega t) (\omega =0.5Hz)');
legend('r', 'y')
xlabel('Time [s]');
ylabel('Amplitude [rad]');

S = feedback(1,openLoop);
PS = plant*S;
disp(norm(PS, inf));
w = logspace(-1,3,1000);
[mag, ~] = bode(PS,w);
mag =squeeze(mag);
figure;
semilogx(w,mag,'Color', [255/255 140/255 105/255]);
hold on;
plot(w, 1/3*ones(size(w)), 'Color', [204/255 51/255 51/255]);
grid on;
title('Magnitude of PS at K = 125 [Absolute Units]');
xlabel('Frequency [rad/s]');
ylabel('Magnitude [abs]');
legend('PS', 'Max Gain 1/3');
hold off;

figure;
impulse(PS);
grid on;
title('Disturbance Impulse Response');
xlabel('Time [s]');
ylabel('Output \phi(t) [rad]');

figure;
step(PS);
grid on;
title('Disturbance Step Response');
xlabel('Time [s]');
ylabel('Output \phi(t) [rad]');