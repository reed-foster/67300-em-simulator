clear all
close all
%clc

%parameters 
p.N = 1000; % discretization
p.eps_0 =  8.85e-12; % F/m
p.mu_0 = 1.26e-6; % N/A^2
p.omega_chi_1 = 2*pi*2.9e15; % Hz
p.chi_1 = 4;
p.chi_2 = 41.7e-12; % m/V %Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
p.delta_chi_1 = 2*pi*1.2e13; % Hz
p.delta_x = 5e-9; % m

% nodal variables
E = zeros(p.N,1); % electric field 
E_dt = zeros(p.N,1); % first derivative of E
H = zeros(p.N,1); % magnetic field 
P = zeros(p.N,1); % polarization vector
P_dt = zeros(p.N,1); % first derivative of P

% source 
J = zeros(p.N,1);

% simulation setup and initial conditions
M = diag([ones(2*p.N, 1); zeros(p.N, 1); ones(2*p.N, 1)]);

X0 = [1e-9*gaussian_start(1,100,400,500,p.N)'; 1e-9*gaussian_start(1,100,400,500,p.N)'; zeros(p.N,1); zeros(p.N,1); zeros(p.N,1)]; % V/m or A/m
Xp0 = OneDimProblem(X0,J,p);
X0(2*p.N+1:3*p.N) = Xp0(p.N+1:2*p.N);
X0(3*p.N+1:4*p.N) = p.eps
X0(4*p.N+1:5*p.N) = Xp0(3*p.N+1:4*p.N);
%X0([1 N N+1 2*N+1]) = 0;
%X0 = [zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1) zeros(N,1)]; % V/m
tspan = linspace(0,1e-12,200); % s

%implicitDAE = @(Xp) M*Xp - M*OneDimProblem(X0,J,p);
%Xp = fsolve(implicitDAE, X0);

options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-4,'AbsTol',1e-8);

[t,X] = ode15s(@(t,X) OneDimProblem(X,J,p), tspan, X0, options);

H_t = X(:,1:p.N);
E_t = X(:,p.N+1:2*p.N);

x = linspace(0,(p.N-1)*delta_x,p.N);
for i=1:10
  figure(i);
  plot(x*1e6, E_t(20*(i-1)+1,:), '-o');
  hold on;
  plot(x*1e6, H_t(20*(i-1)+1,:), '-o');
  
  legend("E_z(x,0)", "H_y(x,0)");
  xlabel("x [um]");
  ylabel("field [A/m or V/m]");
  title(sprintf("t = %d [fs]", tspan(i)));
end
