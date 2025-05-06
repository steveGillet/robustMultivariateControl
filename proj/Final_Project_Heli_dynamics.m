% MCEN 6228 Final Project
%
% The following rigid body model of a rotary-wing MAV has 11 states, 3 inputs, 
% and 3 outputs for controller synthesis. 
%
% The states (in order) are forward speed (u), lateral velocity (v), roll rate (p), pitch 
% rate (q), roll attitude (phi), pitch attitude (theta), a_s and b_s (rotor flap states),
% heave (w), yaw rate (r), yaw rate_feedback (r_fb). The control inputs are lateral and 
% longitudinal cyclic for actuating roll and pitch dynamics (delta_lat, delta_long), 
% and tail rotor (delta_ped) for actuating yaw dynamics. The outputs available 
% for feedback are roll attitude (phi), pitch (theta), and washed out yaw rate 
% (r_fb). 

clear all; 

% Rotorcraft dynamics
A = [-0.1778,zeros(1,4),-9.7807,-9.7807,zeros(1,4);
     0,-0.3104,0,0,9.7807,0,0,9.7807,zeros(1,3);
     -0.3326,-0.5353,zeros(1,4),75.7640,343.86,zeros(1,3);
     0.1903,-0.294,zeros(1,4),172.62,-59.958,zeros(1,3);
     0,0,1,zeros(1,8);zeros(1,3),1,zeros(1,7);
     zeros(1,3),-1,0,0,-8.1222,4.6535,zeros(1,3);
     0,0,-1,zeros(1,3),-0.0921,-8.1222,zeros(1,3);
     zeros(1,6),17.168,7.1018,-0.6821,-0.1070,0;
     0,0,-0.2834,zeros(1,5),-0.1446,-5.5561,-36.674;
     zeros(1,9),2.7492,-11.1120];
 
B = [zeros(6,3);
     0.0632,3.339,0;            % delta_lat (roll motions)
     3.1739,0.2216,0;           % delta_long (pitch motions)
     zeros(1,3);
     0,0,-74.364;               % delta_ped (yaw motions)
     zeros(1,3)];

C = [zeros(1,4),1,zeros(1,6);   % phi
     zeros(1,5),1,zeros(1,5);   % theta
     zeros(1,9),0,1];           % r_fb

D = zeros(size(C,1),size(B,2));

G = ss(A,B,C,D);

