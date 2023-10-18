clear
close all

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

lambda0 = 1.55e-6; % m
eta_0 = 1/sqrt(p.eps_0/p.mu_0);
X0 = [gaussian_start(1/eta_0,2*lambda0,4*lambda0,lambda0,p.N,p.delta_x,0)'; -gaussian_start(1,2*lambda0,4*lambda0,lambda0,p.N,p.delta_x,0)']; % V/m and A/m


tspan = linspace(0,100e-15,331); % s

t_start = 0;
t_stop = 100e-15;
timestep = 10e-18;

[X,t] = ForwardEuler(@evalf_freespace,X0,J,p,t_start,t_stop,timestep,0);

tspan = linspace(t_start, t_stop, size(X,2));
H = X(1:p.N,:)';
E = X(p.N+1:2*p.N,:)';

figure(1);

video = VideoWriter('evalf_freespace_test_Euler.avi'); %Create a video object
open(video); % Open video source - restricts the use of video for your program

for i=1:size(tspan,2)
  yyaxis left;
  plot(x*1e6, E(i,:), '-o');
  ylabel("field [V/m]");
  ylim([-2 2]);
  yyaxis right;
  plot(x*1e6, 1e3*H(i,:), '-o');
  ylim([-10 10]);
  ylabel("field [mA/m]");
  
  legend("E_z(x,0)", "H_y(x,0)");
  xlabel("x [um]");
  title(sprintf("t = %d [fs]", tspan(i)*1e15));
  drawnow;
  vidFrame = getframe(gcf);
  if i < size(tspan,2)
    clf;
  end
  writeVideo(video,vidFrame);
end

close(video);


