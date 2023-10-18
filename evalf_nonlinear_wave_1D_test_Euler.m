close all

params = nonlinear_params();

%parameters 
params.N = 150; % discretization
params.eps_0 =  8.85e-12; % F/m
params.mu_0 = 1.26e-6; % N/A^2
params.delta_x = 100e-9; % m

% nodal variables
E = zeros(params.N,1); % electric field 
H = zeros(params.N,1); % magnetic field 

% source 
J = zeros(params.N,1);

E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz, 1), params.N);
dtP0 = zeros(size(params.Lorentz, 1), params.N);

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.delta_x,params.N); % m


%X0 = [sqrt(params.eps_0/params.mu_0)*gaussian_start(1,100,400,100,params.N)'; -gaussian_start(1,100,400,100,params.N)']; % V/m and A/m
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

tspan = linspace(0,100e-15,331); % s
t_start = 0;
t_stop = 1e-6;
timestop = 10e-9;

%[t,X] = ode45(@(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params), tspan, X0, options);

[X,t] = ForwardEuler(@nonlinear_f, X0, dtJ, p, t_start,t_stop,timestep,0);

H_t = X(:,1:params.N);
E_t = X(:,params.N+1:2*params.N);

figure(1);

video = VideoWriter('evalf_nonlinear_wave_1D_test_Euler.avi'); %Create a video object
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


