clear all
close all

params = nonlinear_params();

% source 
J = zeros(params.N,1);
dtJ = zeros(params.N,1);

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.dz,params.N);

eta = 1/sqrt(params.eps_0*params.Lorentz(1,1)/params.mu_0);
E0 = gaussian_start(1/eta, 1e-6, 5e-6, 1.55e-6, params.N, params.dz);
dtE0 = zeros(params.N, 1);
P0 = params.eps_0*params.Lorentz(:,1)*E0;
dtP0 = zeros(params.N, 1);

X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

tspan = linspace(0,100e-18,1230); % s

options = odeset('RelTol',1e-3,'AbsTol',1e-6);

[t,X] = ode45(@(t,X) nonlinear_f(X,dtJ,params), tspan, X0, options);


figure(1);

video = VideoWriter('nonlinear_f_test.avi'); %Create a video object
open(video); % Open video source - restricts the use of video for your program

for i=1:size(tspan,2)
  [E, dtE, P, dtP] = nonlinear_split_X(X(i,:),params);
  plot(x*1e6, E*1e3, '-o');
  ylabel("field [mV/m]");
  ylim([-10 10]);
  
  legend("E_x(z,t)");
  xlabel("x [um]");
  title(sprintf("t = %0.3f [as]", tspan(i)*1e18));
  drawnow;
  vidFrame = getframe(gcf);
  if i < size(tspan,2)
    clf;
  end
  writeVideo(video,vidFrame);
end

close(video);
