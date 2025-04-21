% F-16 lateral regulator
% JSH 4/10/25

clear all global;

% Plant dynamics
A =[-0.322 0.0640 0.0364 -0.9917 0.0003 0.0008 0;   % x1 = beta
    0 0 1 0.0037 0 0 0;                             % x2 = phi
    -30.6492 0 -3.6784 0.6646 -0.7333 0.1315 0;     % x3 = p
    8.5396 0 -0.0254 -0.4764 -0.0319 -0.0620 0;     % x4 = r
    0 0 0 0 -20.2 0 0;                              % x5 = delta_a
    0 0 0 0 0 -20.2 0;                              % x6 = delta_r
    0 0 0 57.2958 0 0 -1];                          % x7 = x_w

B = [0 0; 0 0; 0 0; 0 0; 20.2 0; 0 20.2; 0 0];  % u1 = delta_a (aileron)
                                                % u2 = delta_r (rudder)

C = [0 0 0 57.2958 0 0 -1;  % y1 = r_w (washed out yaw rate)
     0 0 57.2958 0 0 0 0;   % y2 = p
     57.2958 0 0 0 0 0 0;   % y3 = beta
     0 57.2958 0 0 0 0 0];  % y4 = phi

D = zeros(4,2);

% Static feedback controller u = -Ky
K = [-0.56 -0.44 -0.11 -0.35; -1.19 -0.21 -0.44 0.26];

% Plant dynamics
cA =[0 0 0 0 0 0 0;   
     0 0 0 0 0 0 0;                             
     0 0 0 0 0 0 0;     
     0 0 0 0 0 0 0;     
     0 0 0 0 0 0 0;                              
     0 0 0 0 0 0 0;                              
     0 0 0 0 0 0 0];                          

cB = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
                                            

cC = [0 0 0 0 0 0 0;  
      0 0 0 0 0 0 0;];  

cD = K;

plantSS = ss(A,B,C,D);
controllerSS = ss(cA,cB,cC,cD);

Ti = controllerSS*plantSS*inv(eye(size(plantSS)) + controllerSS*plantSS);
So = inv(eye(size(plantSS)) + plantSS*controllerSS);

s = tf('s');
Wp = (s/10 + 3)/(s+3*4);
Wo = (0.02*s+0.05)/(0.02*s/0.4+1);

N = [-Wo + Ti; -K*So; Wp*G*So; Wp*So];
omega = logspae(-3,3);
muPlotOptions = bodeoptions;
muPlotOptions.PhaseVisible = 'off';
muPlotOptions.XLim = [1e-3 1e3];
muPlotOptions.MagUnits = 'abs';

woTiFRD = frd(-Wo + Ti, omega);
BlkStruct = [2 2];
[muWoTiFull] = mussv(woTiFRD, BlkStruct);
figure;
bodeplot(muWoTiFull, muPlotOptions);