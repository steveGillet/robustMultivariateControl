% A = [0 1; -1 -2];
% B = [0; 1];
% C = [1 0; 0 2];
% D = [0];

A = [-1.5 1 0.1; -4 -1 0; 1 0 0];
B = [-0.2; -1.8; 0];
C = [1 0 0; 0 1 0; 0 0 1];
D = [1; 0; 0];

nStates = size(A,1);
nOutputs = size(C,1);
nInputs = size(B,1);

setlmis([]);

P = lmivar(1, [nStates 1]);

epsilon = 1e-5;

lmiterm([-1 1 1 P], 1, 1);

lmiterm([-2 1 1 0], epsilon);
lmiterm([-2 2 2 0], 1);

lmiterm([3 1 1 P], A', 1, 's');
lmiterm([3 1 1 0], C'*C);
lmiterm([3 1 2 P], 1, B);
lmiterm([3 1 2 0], C'*D);
lmiterm([3 2 2 0], D'*D);
lmiterm([-3 1 1 0], epsilon);
lmiterm([-3 2 2 0], 1);

lmisys = getlmis;

[aa, xopt] = gevp(lmisys, 1);

disp("Gamma: ")
disp(sqrt(aa));

disp("'hinfnorm': ")
sys = ss(A,B,C,D);
disp(hinfnorm(sys));