clear 
close all
clc
%{ The following program will simulate the response of local 
% behavior of the vertical and horiztonal aircraft position around the 
% equilibrium point with a State-Feedback Controller & an Observer
% State Feedback Controller.
%
% State-Feedback Controller 
% - Compute State-Feedback Gain K with State-Feedback desired poles
% - Convert Linearized Controllable Canonical Matrices with Gain to 
%   closed loop state-space form & transfer functions.  
% - Plot bode, impulse, step, and high/low frequency sinusoidal 
%   responses for each transfer function.
%
% Observer State-Feedback Controller 
% - Compute Observer Gain L with Observable desired poles & state-feedback
% Gain with the same controllable desired poles in the State-Feedback
% Controller
% - Convert Linearized Jacobian Matrices with the State-Feedback and 
%   Observer Gains to closed loop state-space form & transfer functions.  
% - Plot bode, impulse, step, and high/low frequency sinusoidal 
%   responses for each transfer function.
%}


%State-Feedback Controller
% Controllable Canonical State-Space
A_cntrl = [0 1 0 0 0 0;
0 0 1 0 0 0;
0 0 0 1 0 0;
0 0 0 0 1 0;
0 0 0 0 0 1;
0 0 0 0 -0.0001563 -0.025];

B_cntrl = [0  0;
0 0;
0 0;
0 0;
0 0;
1 1;];

%1st row is X State-Space & 2nd row is Y State-Space
C_cntrl= [0 -0.6454 -51.63 0.003125 0.25 0; 0 0 0 0.003125 0.25 0];
D_cntrl = [0 0; 0 0];

% Compute State-Feedback Gain
DesPoles_cntrl = [-0.25 + 1j, -0.25 - 1j, -5, -10, -25, -50];
K = place(A_cntrl,B_cntrl,DesPoles_cntrl); 
format shortG
display(K)

% Compute closed loop State-Space & local X & Y transfer functions 
Acl = A_cntrl - B_cntrl*K;
Ccl = C_cntrl  - D_cntrl*K;
[Xnum, Xden] = ss2tf(Acl, B_cntrl,Ccl, D_cntrl,1);
[Ynum, Yden] = ss2tf(Acl, B_cntrl,Ccl, D_cntrl,2);
X_Tf_ctrl = tf(Xnum(1,:), Xden);
Y_Tf_ctrl = tf(Ynum(2,:), Yden);

% Verify eigenvalues of closed-loop system
disp('eigvalues of Close-Loop, A-BK')
disp(eig(Acl))

disp('Local X Transfer Function:')
display(X_Tf_ctrl)

disp('Local Y Transfer Function:')
display(Y_Tf_ctrl)

% Impulse & Step response of local X
figure(1)
subplot(2,1,1)
impulse(X_Tf_ctrl)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Impulse State-Feedback of X local behavior')
grid on

subplot(2,1,2)
step(X_Tf_ctrl)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Step State-Feedback of X local behavior')
grid on

% Bode Plot of local X
figure(2)
subplot(3,1,1)
bode(X_Tf_ctrl)
grid on

% Sinuoidal input for local X 
t_x = linspace(0, 0.1, 100);  % Time Vector for local X

% Low frequency at 0.01 
omega_lw_x = 0.01;
u_X_low = sin(omega_lw_x*t_x); % Forcing Function for local X
XLowfreq = lsim(X_Tf_ctrl, u_X_low, t_x);

% High frequency at 1000 
omega_hi_x = 1000;
u_X_hi = sin(omega_hi_x*t_x); % Forcing Function for local X
XHifreq = lsim(X_Tf_ctrl, u_X_hi , t_x);

% Plot Low Frequency Sinusoidal response for local X
subplot(3,1,2)
plot(t_x,XLowfreq)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Low Freq. State-Feedback of X local behavior')
grid on

subplot(3,1,3)
plot(t_x,XHifreq)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('High Freq. State-Feedback of X local behavior')
grid on


% Local Y Simulation Response 
% Impulse response of local Y
figure(3)
subplot(2,1,1)
impulse(Y_Tf_ctrl)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Impulse State-Feedback of Y local behavior')
grid on

% Step response of local Y
subplot(2,1,2)
step(Y_Tf_ctrl)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Step State-Feedback of Y local behavior')
grid on

% Bode Plot of local Y
figure(4)
subplot(3,1,1)
bode(Y_Tf_ctrl)
grid on

% Sinuoidal input for local Y
t_y = linspace(0, 0.1, 100);  % Time Vector for local Y

% Low frequency at 0.01 
omega_lw_y = 0.01;
u_Y_low = sin(omega_lw_y*t_y); % Forcing Function for local Y
YLowfreq = lsim(Y_Tf_ctrl, u_Y_low, t_y);

% High frequency at 1000 
omega_hi_y = 1000;
u_Y_hi = sin(omega_hi_y*t_y); % Forcing Function for local Y
YHifreq = lsim(Y_Tf_ctrl, u_Y_hi , t_y);

% Plot Low Frequency Sinusoidal response for local Y
subplot(3,1,2)
plot(t_y,YLowfreq)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Low Freq. State-Feedback of Y local behavior')
grid on

% Plot High Frequency Sinusoidal response for local Y
subplot(3,1,3)
plot(t_y,YHifreq)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('High Freq. State-Feedback of Y local behavior')
grid on

%%


%Observer State-Feedback Controller

% Model parameters
J = 0.0475; %kg m^2
m = 4; %kg
r = 0.25; %m
g = 9.81; % m/s^2
c = 0.05; %Ns/m

%Jacobian Matrices 
A = [0 0 0 1 0 0; 
    0 0 0 0 1 0;
    0 0 0 0 0 1; 
    0 0 -g -c/m 0 0;
    0 0 0 0 -c/m 0;
    0 0 0 0 0 0];

B = [0 0;
    0 0;
    0 0;
    1/m 0;
    0 1/m;
    r/J 0];
C = [1 0 0 0 0 0;
    0 1 0 0 0 0];
D = [0 0; 
    0 0];

%Note A,C,D are the same for X and Y transfer functions

% Compute Observable Gain
DesPoles_obs = [-300 + 10j, -300 - 10j, -1500, -2000, -3000, -5000];
L = place(A',C',DesPoles_obs)';

%Note use the same controllable desired poles
K_obs = place (A,B, DesPoles_cntrl);

%Convert to close loop state space & transfer function
Acl_obs = [A -B*K_obs ;L*C (A - L*C- B*K_obs)];
display(L)
display(K_obs)
Bcl_obs = [B; B];
Ccl_obs = [C zeros(2,6)];
Dcl_obs = D;

[Xnum_obs, Xden_obs] = ss2tf(Acl_obs,Bcl_obs, Ccl_obs, Dcl_obs,1);
[Ynum_obs, Yden_obs] = ss2tf(Acl_obs,Bcl_obs, Ccl_obs, Dcl_obs,2);
X_Tf_obs = tf(Xnum_obs(1,:), Xden_obs);
Y_Tf_obs = tf(Ynum_obs(2,:), Yden_obs);


% Verify Closed Loop Poles
disp('eigvalues of Close-Loop, Observer-State Feedback')
disp(eig(Acl_obs))


% Impulse & Step response of local X
figure(5)
subplot(2,1,1)
impulse(X_Tf_obs)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Impulse Observer State-Feedback of X local behavior')
grid on

subplot(2,1,2)
step(X_Tf_obs)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Step Observer State-Feedback of X local behavior')
grid on

% Bode Plot of local X
figure(6)
subplot(3,1,1)
bode(X_Tf_obs)
grid on

% Sinuoidal input for local X 
t_x = linspace(0, 0.1, 100);  % Time Vector for local X

% Low frequency at 0.01 
omega_lw_x = 0.01;
u_X_low = sin(omega_lw_x*t_x); % Forcing Function for local X
XLowfreq_obs = lsim(X_Tf_obs, u_X_low, t_x);

% High frequency at 1000 
omega_hi_x = 1000;
u_X_hi = sin(omega_hi_x*t_x); % Forcing Function for local X
XHifreq_obs = lsim(X_Tf_obs, u_X_hi , t_x);

% Plot Low Frequency Sinusoidal response for local X
subplot(3,1,2)
plot(t_x,XLowfreq_obs)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Low Freq. Observer State-Feedback of X local behavior')
grid on

subplot(3,1,3)
plot(t_x,XHifreq_obs)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('High Freq. Observer State-Feedback of X local behavior')
grid on


% Local Y Simulation Response 
% Impulse response of local Y
figure(7)
subplot(2,1,1)
impulse(Y_Tf_obs)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Impulse Observer State-Feedback of Y local behavior')
grid on

% Step response of local Y
subplot(2,1,2)
step(Y_Tf_obs)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Step Observer State-Feedback of Y local behavior')
grid on

% Bode Plot of local Y
figure(8)
subplot(3,1,1)
bode(Y_Tf_obs)
grid on

% Sinuoidal input for local Y
t_y = linspace(0, 0.1, 100);  % Time Vector for local Y

% Low frequency at 0.01 
omega_lw_y = 0.01;
u_Y_low = sin(omega_lw_y*t_y); % Forcing Function for local Y
YLowfreq_obs = lsim(Y_Tf_obs, u_Y_low, t_y);

% High frequency at 1000 
omega_hi_y = 1000;
u_Y_hi = sin(omega_hi_y*t_y); % Forcing Function for local Y
YHifreq_obs = lsim(Y_Tf_obs, u_Y_hi , t_y);

% Plot Low Frequency Sinusoidal response for local Y
subplot(3,1,2)
plot(t_y,YLowfreq_obs)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Low Freq. Observer State-Feedback of Y local behavior')
grid on

% Plot High Frequency Sinusoidal response for local Y
subplot(3,1,3)
plot(t_y,YHifreq_obs)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('High Freq. Observer State-Feedback of Y local behavior')
grid on
