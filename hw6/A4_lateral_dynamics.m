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
