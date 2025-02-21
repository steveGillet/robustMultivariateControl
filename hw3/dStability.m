A = [-1.5 1 0.1; -4 -1 0; 0 1 0];
% A = [0 1; -1 -2];
nStates = size(A, 1);
disp(eig(A));

setlmis([])

P = lmivar(1, [nStates 1]);
alpha = 0.05;

lmiterm([-1 1 1 P],1,1);
lmiterm([2 1 1 P],1,A,'s');
lmiterm([2 1 1 P],2*alpha,1);

lmiDstab = getlmis;

[tmin, pfeas] = feasp(lmiDstab);
P = dec2mat(lmiDstab,pfeas,P);
disp(P);
