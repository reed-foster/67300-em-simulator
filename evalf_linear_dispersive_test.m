clear all
close all
%clc

%parameters 
p.N = 1500; % discretization
p.eps_0 =  8.85e-12; % F/m
p.mu_0 = 1.26e-6; % N/A^2
p.chi_1 = 4;
p.omega_chi_1 = 2*pi*2.9e9; % Hz
p.delta_chi_1 = 2*pi*1.2e6; % Hz
p.delta_x = 100e-9; % m

% nodal variables
E = zeros(p.N,1); % electric field 
H = zeros(p.N,1); % magnetic field 
P = zeros(p.N,1); % polarization vector
P_dt = zeros(p.N,1); % first derivative of P

% source 
J = zeros(p.N,1);

% simulation setup and initial conditions
x = linspace(0,(p.N-1)*p.delta_x,p.N);

M = eye(4*p.N);

X0 = [sqrt(p.eps_0*(1+p.chi_1)/p.mu_0)*gaussian_start(1,100,400,100,p.N)'; -gaussian_start(1,100,400,100,p.N)'; zeros(p.N,1); zeros(p.N,1)]; % V/m or A/m
X0(2*p.N+1:3*p.N) = p.eps_0*p.chi_1*X0(p.N+1:2*p.N);

Xp0 = evalf_linear_dispersive(X0,J,p);

%figure;
%subplot(2,1,1);
%yyaxis left;
%plot(x*1e6, X0(1:p.N));
%yyaxis right;
%plot(x*1e6, X0(p.N+1:2*p.N));
%legend("H_y(x,0)", "E_z(x,0)");
%subplot(2,1,2);
%yyaxis left;
%plot(x*1e6, Xp0(1:p.N));
%yyaxis right;
%plot(x*1e6, Xp0(p.N+1:2*p.N));
%legend("d_tH_y(x,0)", "d_tE_z(x,0)");
%
%figure;
%yyaxis left;
%plot(x*1e6, Xp0(2*p.N+1:3*p.N));
%yyaxis right;
%plot(x*1e6, Xp0(3*p.N+1:4*p.N));
%legend("d_tP_z(x,0)", "d_tDP_z(x,0)");

tspan = linspace(0,100e-15,1001); % s

options = odeset('Mass',M,'MassSingular','no','RelTol',1e-3,'AbsTol',1e-6);

[t,X] = ode45(@(t,X) evalf_linear_dispersive(X,J,p), tspan, X0, options);

H_t = X(:,1:p.N);
E_t = X(:,p.N+1:2*p.N);

figure(1);

video = VideoWriter('evalf_linear_dispersive_test.avi'); %Create a video object
open(video); % Open video source - restricts the use of video for your program
for i=1:size(tspan,2)
  yyaxis left;
  plot(x*1e6, E_t(i,:), '-o');
  ylabel("field [V/m]");
  ylim([-4 4]);
  yyaxis right;
  plot(x*1e6, 1e3*H_t(i,:), '-o');
  ylim([-20 20]);
  ylabel("field [mA/m]");
  
  legend("E_z(x,0)", "H_y(x,0)");
  xlabel("x [um]");
  title(sprintf("t = %d [fs]", tspan(i)*1e15));

  drawnow;
  vidFrame = getframe(gcf);
  clf;
  writeVideo(video,vidFrame);
end

close(video);
