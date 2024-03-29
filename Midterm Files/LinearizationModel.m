clear 
close all
clc
%{ The following program will simulate the local behavior of the vertical 
% and horiztonal aircraft position around the equilibrium point.
%
% - Convert Linearized Jacobian Matrices to state-space form & transfer 
%   functions. 
% - Compute the poles and zeros & plot bode for each transfer function. 
% - Plot the impulse, step, and high/low frequency sinusoidal 
%   responses for each transfer function.
%}

% Model parameters
J = 0.0475; %kg m^2
m = 4; %kg
r = 0.25; %m
g = 9.81; % m/s^2
c = 0.05; %Ns/m

% Linearized Jacobian Matrices evaluated at the equilibrium point
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

% Convert Jacobian Matrices to State-Space Form
stateSpace = ss(A,B,C,D);
% Visualize the State-Space Form
disp('State-Space Form:')
display(stateSpace)

%%
% Local X Simulation Response 
% Convert State-Space Form to local X Transfer Function
[X_num, X_den] = ss2tf(A,B,C,D,1);
X_tf = tf(X_num(1,:),X_den);

% Visualize the local X transfer function
disp('local X Transfer Function:')
display(X_tf)

% Visualize local X Transfer Function Poles
disp('X Transfer Function Poles:')
display(pole(X_tf));

% Visualize local X Transfer Function Zeros
disp('local X Transfer Function Zeros:')
display(zero(X_tf));

% Impulse response of local X
figure(1)
impulse(X_tf)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Impulse response of local behvaior of X')
grid on

% Step response of local X
figure(2)
step(X_tf)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Step response of X local behvaior')
grid on

% Bode Plot of local X
figure(3)
bode(X_tf)
grid on

% Sinuoidal input for local X 
t_x = linspace(0, 1000, 1000);  % Time Vector for local X

% Low frequency at 0.01 
omega_lw_x = 0.01;
u_X_low = sin(omega_lw_x*t_x); % Forcing Function for local X
XLowfreq = lsim(X_tf, u_X_low, t_x);

% High frequency at 1000 
omega_hi_x = 1000;
u_X_hi = sin(omega_hi_x*t_x); % Forcing Function for local X
XHifreq = lsim(X_tf, u_X_hi , t_x);

% Plot Low Frequency Sinusoidal response for local X
figure(4)
plot(t_x,XLowfreq)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('Low Frequency Sinusoidal response of X local behvaior')
grid on

% Plot High Frequency Sinusoidal response for local X
figure(5)
plot(t_x,XHifreq)
xlabel('time (s)')
ylabel('x - x_e position (m)')
title('High Frequency Sinusoidal X local behvaior')
grid on


%%
% Local Y Simulation Response 
% Convert State-Space Form to Y Transfer Function
[Y_num, Y_den] = ss2tf(A,B,C,D,2);
disp('Y Transfer Function:')
Y_tf = tf(Y_num(2,:),Y_den);

% Visualize the Y transfer function
disp('Y Transfer Function:')
display(Y_tf)

% Visualize Y Transfer Function Poles
disp('Y Transfer Function Poles:')
display(pole(Y_tf));

% Visualize Y Transfer Function Zeros
disp('Y Transfer Function Zeros:')
display(zero(Y_tf));

% Impulse response of local Y
figure(6)
impulse(Y_tf)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Impulse response of Y local behvaior')
grid on

% Step response of local Y
figure(7)
step(Y_tf)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Step response of Y local behvaior')
grid on

% Bode Plot of local Y
figure(8)
bode(Y_tf)
grid on

% Sinuoidal input for local Y
t_y = linspace(0, 1000, 1000);  % Time Vector for local Y

% Low frequency at 0.01 
omega_lw_y = 0.01;
u_Y_low = sin(omega_lw_y*t_y); % Forcing Function for local Y
YLowfreq = lsim(Y_tf, u_Y_low, t_y);

% High frequency at 1000 
omega_hi_y = 1000;
u_Y_hi = sin(omega_hi_y*t_y); % Forcing Function for local Y
YHifreq = lsim(Y_tf, u_Y_hi , t_y);

% Plot Low Frequency Sinusoidal response for local Y
figure(9)
plot(t_y,YLowfreq)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('Low Frequency Sinusoidal response of Y local behvaior')
grid on

% Plot High Frequency Sinusoidal response for local Y
figure(10)
plot(t_y,YHifreq)
xlabel('time (s)')
ylabel('y - y_e position (m)')
title('High Frequency Sinusoidal of Y local behvaior')
grid on
