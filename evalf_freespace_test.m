clear all
close all
%clc

%parameters 
p.N = 1500; % discretization
p.eps_0 =  8.85e-12; % F/m
p.mu_0 = 1.26e-6; % N/A^2
p.delta_x = 100e-9; % m

% nodal variables
E = zeros(p.N,1); % electric field 
H = zeros(p.N,1); % magnetic field 

% source 
J = zeros(p.N,1);

% simulation setup and initial conditions
x = linspace(0,(p.N-1)*p.delta_x,p.N);

M = kron([1 0; 0 1], eye(p.N));

X0 = [sqrt(p.eps_0/p.mu_0)*gaussian_start(1,100,400,100,p.N)'; -gaussian_start(1,100,400,100,p.N)']; % V/m and A/m

Xp0 = evalf_freespace(X0,J,p);

figure(2);
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

tspan = linspace(0,100e-15,331); % s

options = odeset('Mass',M,'MassSingular','no','RelTol',1e-3,'AbsTol',1e-6);

[t,X] = ode45(@(t,X) evalf_freespace(X,J,p), tspan, X0, options);

H_t = X(:,1:p.N);
E_t = X(:,p.N+1:2*p.N);

figure(1);

video = VideoWriter('evalf_freespace_test.avi'); %Create a video object
open(video); % Open video source - restricts the use of video for your program

for i=1:size(tspan,2)
  yyaxis left;
  plot(x*1e6, E_t(i,:), '-o');
  ylabel("field [V/m]");
  ylim([-2 2]);
  yyaxis right;
  plot(x*1e6, 1e3*H_t(i,:), '-o');
  ylim([-10 10]);
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
