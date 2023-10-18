close all

params = nonlinear_params();

% nodal variables
E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz, 1), params.N);
dtP0 = zeros(size(params.Lorentz, 1), params.N);
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.dz,params.N); % m

t_start = 0;
t_stop = 80e-15;
dt = 1e-18; % min 1e-18 max 5e18
tsteps = round(t_stop/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s

% source 
dtJ = zeros(params.N,tsteps);
omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
t0_J = 1.5*2*pi/omega_J;
ampl = 1e5; % also try 1e7, 1e6, 1e8
dtJ(round(params.N/2),:) = ricker(ampl/params.dz*omega_J, omega_J, tsteps, dt, t0_J);

[X,t] = ForwardEuler(@nonlinear_f, X0, dtJ, params, t_start, t_stop, dt, 0);


gen_video = true;
plot_jacobian = false;
x = linspace(0,(params.N-1)*params.dz,params.N); % m
zoom_width = 10e-6; % m
frame_dec = 500;
nonlinear_f_plot(x, tspan, X, params, zoom_width, gen_video, plot_jacobian, frame_dec, 'nonlinear_f_test.avi');

