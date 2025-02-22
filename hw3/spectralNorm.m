% A = [2 -1; -1 2];
A = [-1 1 0; 1 -1 1; 0 1 1];
disp("Singular Values: ");
disp(svd(A));

nStates = size(A,1);

setlmis([])

rho = lmivar(1, [1 0]);

lmiterm([-1 1 1 rho],1,1);
lmiterm([-1 1 2 0],A');
lmiterm([-1 2 2 0],1);

spectralLmi = getlmis;

c = mat2dec(spectralLmi,rho);

[copt, xopt] = mincx(spectralLmi, c);
rhoVal = dec2mat(spectralLmi, xopt, 1);
gamma = sqrt(rhoVal);
disp("Gamma: ");
disp(gamma);
