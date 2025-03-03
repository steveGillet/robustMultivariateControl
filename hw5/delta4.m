num = [40 1];
den = conv([1 1 4], [1 6]);
plant = tf(num, den);

bode(plant);

numC = [100000*0.08 10000];
denC = [100000*0.08 1];
C = tf(numC, denC);
openLoop = plant*C;
[~, ~, ~, Wc] = margin(openLoop);
figure;
bode(openLoop);
hold on;
xline(Wc, 'r-', 'Gain Crossover');