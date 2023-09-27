clear all
close all
%clc

%parameters 
p.N = 1500; % discretization
p.eps_0 =  8.85e-12; % F/m
p.mu_0 = 1.26e-6; % N/A^2
p.delta_x = 5e-9; % m

% nodal variables
E = zeros(p.N,1); % electric field 
H = zeros(p.N,1); % magnetic field 

% source 
J = zeros(p.N,1);

% simulation setup and initial conditions
x = linspace(0,(p.N-1)*p.delta_x,p.N);

M = eye(2*p.N);

X0 = [sqrt(p.eps_0/p.mu_0)*gaussian_start(1,100,400,300,p.N)'; -sqrt(p.mu_0/p.eps_0)*gaussian_start(1,100,400,300,p.N)']; % V/m and A/m

Xp0 = evalf_freespace(X0,J,p);

figure(11);
subplot(2,1,1);
yyaxis left;
plot(x*1e6, X0(1:p.N));
yyaxis right;
plot(x*1e6, X0(p.N+1:2*p.N));
legend("H_y(x,0)", "E_z(x,0)");
subplot(2,1,2);
yyaxis left;
plot(x*1e6, Xp0(1:p.N));
yyaxis right;
plot(x*1e6, Xp0(p.N+1:2*p.N));
legend("d_tH_y(x,0)", "d_tE_z(x,0)");

%return;

tspan = linspace(0,1e-12,10000); % s

options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-4,'AbsTol',1e-8);

[t,X] = ode15s(@(t,X) evalf_freespace(X,J,p), tspan, X0, options);

H_t = X(:,1:p.N);
E_t = X(:,p.N+1:2*p.N);

for i=1:10
  figure(i);
  plot(x*1e6, E_t(1000*(i-1)+1,:), '-o');
  hold on;
  plot(x*1e6, H_t(1000*(i-1)+1,:), '-o');
  
  legend("E_z(x,0)", "H_y(x,0)");
  xlabel("x [um]");
  ylabel("field [A/m or V/m]");
  title(sprintf("t = %d [ps]", tspan(1000*(i-1)+1)*1e12));
end
