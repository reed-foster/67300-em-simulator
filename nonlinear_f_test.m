clear all
close all

params = nonlinear_params();

% source 
params.tsteps = round(params.tf/params.dt);
params.tvec = linspace(0,(params.tsteps-1)*params.dt,params.tsteps); % s

dtJ = @(t) ricker(t, params.ampl_J, params.omega_J, params.t0_J);

% simulation setup and initial conditions
E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz,1), params.N);
dtP0 = zeros(size(params.Lorentz,1), params.N);
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

options = odeset('RelTol',1e-3,'AbsTol',1e-12);
eval_f = @(t,X) nonlinear_f(X,nonlinear_u(t,params),params);
[t,X] = ode45(eval_f, params.tvec, X0, options);

gen_video = true;
plot_jacobian = false;
xvec = linspace(0,(params.N-1)*params.dz,params.N); % m
zoom_width = 10e-6; % m
frame_dec = 10;
nonlinear_f_plot(xvec, params.tvec, X', params, zoom_width, gen_video, plot_jacobian, frame_dec, 'nonlinear_f_test.avi');
