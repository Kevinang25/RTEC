% ********** Q3 CODE **********
clear all
close all
clc
% System parameters
wo = 1;
am = 1;
theta = 0;
d = 1;
deltat = 0.05;
Nsim = 600;
t = 0:deltat:(Nsim-1)*deltat;

% --- Q3.1. State Space Model ---
A = [0 1; -wo^2 0];
C = [1 0];
N = [0 wo^2];

% --- Q3.2. Observer Design---
OL = [C; C*A]; % Observability matrix
Pob = conv([1 5], [1 5]);
Kob = T2place(A',C',Pob)';

% --- Q3.3. Simulation ---
% Generate signals
y0 = am*sin(wo*t+theta); % Signal with no noise
epsilon = 0.2*randn(1,Nsim);
y = y0 + epsilon; % Signal with noise
Xhat = [0;0]; % Assume initial condition is 0.

% Simulation
for k = 1:Nsim
    Xhat = Xhat+deltat*(A*Xhat+N*d+Kob*(y(k)-C*Xhat));
    Xhat1(k) = Xhat(1,1);
    Xhat2(k) = Xhat(2,1);
end

figure(1)
%Xhat1 is the estimated signal.
subplot(2,1,1)
plot(t,y,'r',t,Xhat1,'b')
title('Estimated vs. Noisy Signal')
legend('Noisy: y', 'Estimated: $$\hat{x}$$','Interpreter','Latex')
xlabel('Time (s)')
ylabel('Amplitude')
axis([0 30 -1.5 1.5]);
grid on

subplot(2,1,2)
plot(t,y0,'m',t,Xhat1,'b')
title('Estimated vs. Noise-Free Signal')
legend('Noise-free: y0', 'Estimated: $$\hat{x}$$','Interpreter','Latex')
xlabel('Time (s)')
ylabel('Amplitude')
axis([0 30 -1.5 1.5]);
grid on