clear 
close all
clc
%{ The following program will determine the stability of the 
% Linearized System with Eigenvalues and the Lyapunov equation. 
% - The linearized Jacobian Matrices are converted to Jordan Form, 
% Transfer functions, and Discrete State Space
% - The controllability/observability are checked for this system. 
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
disp('Linearized State-Space Form:')
display(stateSpace)



%%
% Stability Analysis
% Eigenvalue Analysis
% (A-lambda*I) = 0
EigValues = eig(A);
disp(' ')
disp('Eigenvalues of A')
disp(EigValues)

%Lyapunvo Analysis
% A^T*P + P*A= Q
Q = eye(6);
disp('Q Matrix')
disp(Q)
try
    P = lyap(A',Q);
    disp('P Matrix')
    disp(P)
    
    disp('Display P')
    disp(P)
    disp('Check if P Transpose = P')
    disp(P')
    disp("these matrices are the same")
    
    disp(' ')
    %check positive definite matrices
    x = [1; -1; 1];
    Vy = x'*P*x;
    derV = x'*(A'*P +P*A)*x;
    disp('Check if lyapunov function is > 0')
    disp(Vy)
    disp('Check if derivative of lyapunov function is < 0')
    disp(derVy)
    disp('Check eigenvalues of P are > 0')
    disp(eig(P))
    disp('Check determinant of P is > 0')
    disp(det(P))
    disp('Check P_11 > 0')
    disp(P(1,1))
    if (P(1,1) > 0 && Vy > 0 && det(P) > 0 && P(1,1) > 0 && derVy < 0)
        for i = eig(P)
            if i < 0
                break
            end
        end
       disp('This is asymptotically stable system because P is positive-definite ')
    else
      disp('This is asymptotically stable system because P is not positive-definite ')
    end
catch error 
    % if error with lyapunov does not exist
    disp(error.message)
end 
%%
% Jordan Form
% J = V^-1*A*V
[V_g, J_A] = jordan(A);
disp(' ')
disp('Generalized eigenvector Matrix, V_g')
format shortG
disp(round(V_g))

J_B = inv(V_g)*B;
J_C = C* V_g;
J_D = D;

JordForm_ss = ss(J_A,J_B,J_C,J_D);
disp('Jordan form State-Space')
display(JordForm_ss)
%%
%Transfer functions

% Convert Jacobian State-Space Form to local X Transfer Function
[X_num, X_den] = ss2tf(A,B,C,D,1);
X_tf = tf(X_num(1,:),X_den);
disp('Local X Transfer Function:')
display(X_tf)

% Visualize local X Transfer Function Poles
disp('X Transfer Function Poles:')
display(pole(X_tf));

% Visualize local X Transfer Function Zeros
disp('local X Transfer Function Zeros:')
display(zero(X_tf));


% Convert State-Space Form to Y Transfer Function
[Y_num, Y_den] = ss2tf(A,B,C,D,2);
Y_tf = tf(Y_num(2,:),Y_den);
disp('Local Y Transfer Function:')
display(Y_tf)

% Visualize Y Transfer Function Poles
disp('Y Transfer Function Poles:')
disp(pole(Y_tf));

% Visualize Y Transfer Function Zeros
disp('Y Transfer Function Zeros:')
disp(zero(Y_tf));


%%
%Discrete Jacobian State Space
dt = 0.1;
SS_Discr = c2d(stateSpace,dt,'zoh');
disp('Discrete Jacobian State-Space')
display(SS_Discr)

%%
%Controllability & Observability:
if rank(ctrb(A,B)) == size(A,2)
    disp('The system is full rank and controllable')
end
if rank(obsv(A,C)) == size(C,2)
    disp('The system is full rank and observable')
end

disp(['This system is both observable & controllable so this system will be used' ...
    'for the Observer State-Feedback Controller'])



