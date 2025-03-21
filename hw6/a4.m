%%% GIVEN DYNAMICS CODE, CREDIT DR. HUMBERT %%%

% A-4 Lateral dynamics

clear all global; 

% Lateral stability derivatives and aircraft parameters
Y_b = 43.72;      % ft/s^2
Y_r = 0;
N_b = 4.395 ;     % 1/s^2
N_r = -0.744;     % 1/s
Y_dr = 12.17;     % ft/s^2
N_dr = -4.495;    % 1/s^2
N_da = -0.217;    % 1/s^2
u_0 = 176;        % ft/s

% Aircraft model 
As = [ Y_b/u_0 -(1-Y_r/u_0) ;  N_b   N_r ];  % x1 = beta, x2 = r
Bs = [ 0 Y_dr/u_0  ; N_da N_dr  ];           % u1 = delta_a, u2 = delta_r
Cs = eye(2);                                 % y1 = beta, y2 = r
Ds = [ 0 0 ; 0  0];
G = ss(As,Bs,Cs,Ds);

%%% CUSTOM CODE %%%
A=0.05;
omegaB=1;
M=2;
s = tf('s');
Wp = (s/M+omegaB)/(s+omegaB*A)*eye(2);
W2=[];
W3=[];

[K,CL,GAM,INFO] = mixsyn(G,Wp,W2,W3);
disp(GAM);

L = G*K;
S = inv(eye(2)+L);
T = L*S;

figure(1); clf;
sigma(G, 'r', K, 'b', G*K, 'm');
grid on; legend('\sigma(G)','\sigma(K)','\sigma(G*K)'); title('Open Loop Singular Values');

figure(2); clf;
sigma(S, 'm', inv(Wp), 'r--');
grid on; legend('\sigma(So)','\sigma(Wp)'); title('\sigma(So) vs. |1/Wp(jw)|');

figure(3); clf;
step(T); grid on;
title('Step Input Response');

figure(4); clf;
step(S); grid on;
title('Step Disturbance Response');
