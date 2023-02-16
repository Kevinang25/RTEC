% ********** Q2 CODE **********
clear all
close all
clc
% System parameters
deltat = 0.03;
Nc = 6;
Np = 30;
Nsim = 80;
rw = 0.1;
R_bar = rw*eye(Nc,Nc);

% Converting continuous time system to discrete time system
num = 6;
den = [1 0.1 6];
G = tf(num, den);
H = c2d(G, deltat);
% Discrete time system in state space model
[Ap,Bp,Cp,Dp] = tf2ss(num,den);
[Am,Bm,Cm,Dm] = c2dm(Ap,Bp,Cp,Dp,deltat);
[Phi_Phi,Phi_F,Phi_R,A,B,C]=mpcgain(Am,Bm,Cm,Nc,Np);

% --- Q2.1. Eigenvalues of discrete time system ---
op_poles = eig(Am);
figure(1)
plot(op_poles,'x');
title('Open-Loop Poles Location');
zgrid;
% Check for stability
abs(op_poles(1));
abs(op_poles(2));

% --- Q2.2. Observer Design ---
OL = [C; C*A; C*A*A];
rank(OL);
Pob = [0.1 0.2 0.3];
Kob = place(A',C',Pob)';

% --- Q2.3. MPC Design ---
Kmpc_temp = inv(Phi_Phi+R_bar)*Phi_F;
Kmpc = Kmpc_temp(1,:);
Acl = A - B*Kmpc;
cl_poles = eig(Acl);
figure(2)
plot(cl_poles,'x');
title('Closed-Loop Poles Location');
zgrid;

% --- Q2.4. Closed-Loop Simulation ---
% Signal definitions
ref = -2;
y_k = zeros(size(Cm,1),1);
x_m = zeros(size(Am,2),1);
u_k = zeros(size(Bm,2),1);
Xhat = zeros(size(A,1),1);
disturbance = 2*[zeros(1,Nsim/2) ones(1,Nsim/2)];

% Closed-loop simulation
for k = 1:Nsim
    
    Phi_Rs = Phi_R*ref;
    Delta_U = inv(Phi_Phi+R_bar)*(Phi_Rs-Phi_F*Xhat); 
    delta_u_k = Delta_U(1); % Receding Horizon, using first sample
    
    % Observer Implementation
    Xhat = A*Xhat + B*delta_u_k + Kob*(y_k - C*Xhat); 

    u_k = u_k + delta_u_k; % add incremental 'delta_u' to 'u_past' 
     
    % Actual plant
    x_m = Am*x_m + Bm*(u_k + disturbance(k)); % Discretised model.
    y_k = Cm*x_m;

    %Storing the plots
    y_plot(k) = y_k;
    r_plot(k) = ref;
    u_plot(k) = u_k;
end

% Plotting the signals
k_plot = 0:1:(Nsim-1);

figure(3)
subplot(2,1,1)
plot(k_plot,y_plot,'-r')
hold on
plot(k_plot,r_plot,'--b')
title('Closed-Loop Response')
legend('Output: y', 'Reference: r')
xlabel('Sample')
ylabel('Amplitude')
axis([0 80 -2.5 0.5])
grid on

subplot(2,1,2)
plot(k_plot,u_plot)
legend('Control: u')
xlabel('Sample')
ylabel('Amplitude')
axis([0 80 -10 4])
grid on

% --- Q2.5. & Q2.6. Constraints for MPC ---
%  -0.6 <= delta_u(1) <= 0.4
% -4 <= u(1) <= 5
M = [1 zeros(1,(Nc-1));-1 zeros(1,(Nc-1));1 zeros(1,(Nc-1));-1 zeros(1,(Nc-1))];

% Signal definitions
y_k = zeros(size(Cm,1),1);
x_m = zeros(size(Am,2),1);
u_k = zeros(size(Bm,2),1);
Xhat = zeros(size(A,1),1);

% Closed-loop simulation
for k = 1:Nsim
    
    % gamma matrix is defined in the loop as it depends on the value of u_k
    gamma = [0.4; 0.6; 5-u_k; 4+u_k];
    
    Phi_Rs = Phi_R*ref;
    Delta_U = QPhild(Phi_Phi + R_bar, Phi_F*Xhat - Phi_Rs, M, gamma);
    delta_u_k = Delta_U(1); % Receding Horizon, using first sample
    
    % Observer Implementation
    Xhat = A*Xhat + B*delta_u_k + Kob*(y_k - C*Xhat); 

    u_k = u_k + delta_u_k; % add incremental 'delta_u' to 'u_past' 
    
    % Actual plant
    x_m = Am*x_m + Bm*(u_k + disturbance(k)); % Discretised model.
    y_k = Cm*x_m;

    %Storing the plots
    y_plot(k) = y_k;
    r_plot(k) = ref;
    u_plot(k) = u_k;
    delta_u_plot(k) = delta_u_k;
end

% Plotting the signals
u_max = 5*(ones(1,Nsim));
u_min = -4*(ones(1,Nsim));
deltau_max = 0.4*(ones(1,Nsim));
deltau_min = -0.6*(ones(1,Nsim));

figure(4)
subplot(3,1,1)
plot(k_plot,y_plot,'-r',k_plot,r_plot,'--b')
title('Closed-Loop Response w/ Constraints')
legend('Output: y', 'Reference: r')
xlabel('Sample')
ylabel('Amplitude')
axis([0 80 -3 0.5])
grid on

subplot(3,1,2)
plot(k_plot,u_plot,k_plot,u_max,'--m',k_plot,u_min,'--m')
legend('Control: u', 'Constraint')
xlabel('Sample')
ylabel('Amplitude')
axis([0 80 -6 7])
grid on

subplot(3,1,3)
plot(k_plot,delta_u_plot,'g',k_plot,deltau_max,'--m',k_plot,deltau_min,'--m')
legend('Inc. Control: \Deltau', 'Constraint')
xlabel('Sample')
ylabel('Amplitude')
axis([0 80 -0.8 0.6])
grid on
