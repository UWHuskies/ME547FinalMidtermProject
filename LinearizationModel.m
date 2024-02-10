clear 
close all
clc

%Model parameters
J = 0.0475; %kg m^2
m = 4; %kg
r = 0.25; %m
g = 9.81; % m/s^2
c = 0.05; %Ns/m

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

stateSpace = ss(A,B,C,D);
display(stateSpace)


t = linspace(0, 100, 100);  % Time Vector
u = sin(t); % Forcing Function

[num1, den1] = ss2tf(A,B,C,D,1);
tf1 = tf(num1(1,:),den1);
display(tf1)
tf1Poles = pole(tf1);
display(tf1Poles)
tf1Zeros = zero(tf1);
display(tf1Zeros)

y1 = lsim(tf1, u, t); 

[num2, den2] = ss2tf(A,B,C,D,2);
tf2 = tf(num2(2,:),den2);
display(tf2)
tf2Poles = pole(tf2);
display(tf2Poles)
tf2Zeros = zero(tf2);
display(tf2Zeros)


y2 = lsim(tf2, u, t); 

figure(1)
impulse(tf1)
figure(2)
impulse(tf2)
figure(3)
step(tf1)
figure(4)
step(tf2)

figure(5)
plot(t,y1)
figure(6)
plot(t,y2)

figure(7)
bode(tf1)

figure(8)
bode(tf2)