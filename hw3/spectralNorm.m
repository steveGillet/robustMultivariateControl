% A = [2 -1; -1 2];
A = [-1 1 0; 1 -1 1; 0 1 1];

nStates = size(A,1);

setlmis([])

rho = lmivar(2, [1 1]);

lmiterm([-1 1 1 rho],1,1);
lmiterm([2 1 1 rho],1,1);
lmiterm([2 1 2 0],A');
lmiterm([2 2 1 0],A);
lmiterm([2 2 2 0],1);

spectralLmi = getlmis;

[copt, xopt] = mincx(spectralLmi, rho);
rhoVal = dec2mat(spectralLmi, xopt, rho);
gamma = sqrt(rhoVal);
