clear
close all
clc

% Call the nonlinear_params function
params = nonlinear_params;

p.N = params.N;
p.eps_0 = params.eps_0;
p.mu_0 = params.mu_0;
p.chi_2 = params.chi_2;
p.chi_3 = params.chi_3;
p.Lorentz = params.Lorentz;
p.dz = params.dz;
p.x_order = params.x_order;
p.delta_x = 100e-9; % m

% nodal variables

E = zeros(p.N,1); % electric field
H = zeros(p.N,1); % magnetic field

% source
dt = p.dz/3e8; %50e-18; % s
tf = 500e-15; % s
tsteps = round(tf/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s
J = zeros(p.N,tsteps);
dtJ = zeros(p.N,tsteps);

omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
dtJ_ricker = 2e8/p.dz*omega_J*(1-((tspan-2*pi/omega_J)*omega_J).^2).*exp(-(tspan-2*pi/omega_J).^2/(2*(1/omega_J)^2));
dtJ(round(p.N/2),:) = dtJ_ricker;

E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz, 1), p.N);
dtP0 = zeros(size(p.Lorentz, 1), p.N);

% simulation setup and initial conditions

x = linspace(0,(p.N-1)*p.delta_x,p.N);
M = kron([1 0; 0 1], eye(p.N));
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

%tspan = linspace(0,100e-15,331); % s
%
% options = odeset('Mass',M,'MassSingular','no','RelTol',1e-3,'AbsTol',1e-6);
%
% [t,X] = ode45(@(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params), tspan, X0, options)

t_start = 0;
t_stop = (tsteps-1)*dt;
%timestep = 10e-9;

%[t,X] = ode45(@(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params), tspan, X0, options);

eval_f = @(X, dtJ, p) nonlinear_f(X, dtJ(:, round(t/dt + 0.5)), params);
[X, t] = ForwardEuler(eval_f, X0, J, params, t_start,t_stop,tsteps,0);

%tspan = linspace(t_start, t_stop, size(X,2));

H = X(1:p.N,:)';
E = X(p.N+1:2*p.N,:)';

figure(1);

video = VideoWriter('evalf_nonlinear_wave_1D_test_Euler_trial.avi'); %Create a video object
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
  if i < size(tspan,2)
    clf;
  end
  writeVideo(video,vidFrame);
end

close(video);
