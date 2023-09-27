clear all
close all
%clc

% nodal variables
N = 500; % discretization
E = zeros(N,1); % electric field 
E_dt = zeros(N,1); % first derivative of E
H = zeros(N,1); % magnetic field 
P = zeros(N,1); % polarization vector
P_dt = zeros(N,1); % first derivative of P

%parameters 
eps_0 =  8.85e-12; % F/m
mu_0 = 1.26e-6; % N/A^2
omega_chi_1 = 2*pi*2.9e15; % Hz
chi_1 = 4;
chi_2 = 41.7e-12; % m/V %Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
delta_chi_1 = 2*pi*1.2e13; % Hz
delta_x = 50e-9; % m

% source 
J = zeros(N,1);

% simulation setup and initial conditions
M = diag([ones(2*N, 1); zeros(N, 1); ones(2*N, 1)]);
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-8);

X0 = [zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1)]; % V/m
tspan = linspace(0,100e-15,50); % s

[t,X] = ode15s(@(t,X) f(t,X,J,N,delta_x,eps_0,mu_0,chi_1,chi_2,omega_chi_1,delta_chi_1),tspan,X0,options);

H_t = X(:,1:N);
E_t = X(:,N+1:2*N);

x = linspace(0,(N-1)*delta_x,N);
plot(x*1e6, E_t(10,:), '-o');
hold on;
plot(x*1e6, H_t(10,:), '-o');

legend("E_z(x,0)", "H_y(x,0)");
xlabel("x [um]");
ylabel("field [A/m or V/m]");


function out = f(t,X,J,N,delta_x,eps_0,mu_0,chi_1,chi_2,omega_chi_1,delta_chi_1)
  out = OneDimProblem(X(1:N),X(N+1:2*N),X(2*N+1:3*N),X(3*N+1:4*N),X(4*N+1:5*N),J,N,delta_x,eps_0,mu_0,chi_1,chi_2,omega_chi_1,delta_chi_1);
end
