clear all
close all

params = nonlinear_params();
params.N = 2000;

%figure(3);
%plot_permittivity(params);
%return;

% source 
dt = params.dz/3e8; %50e-18; % s
tf = 300e-15; % s
tsteps = round(tf/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s
J = zeros(params.N,tsteps);
dtJ = zeros(params.N,tsteps);

omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
t0_J = 2*2*pi/omega_J;
dtJ(round(params.N/2),:) = ricker(5e7/params.dz*omega_J, omega_J, tsteps, dt, t0_J);
%figure(4);
%plot(tspan, dtJ(round(params.N/2),:), '-o');
%return;

% simulation setup and initial conditions
E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz,1), params.N);
dtP0 = zeros(size(params.Lorentz,1), params.N);
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

options = odeset('RelTol',1e-3,'AbsTol',1e-12);
eval_f = @(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params);
[t,X] = ode45(eval_f, tspan, X0, options);

gen_video = true;
plot_jacobian = false;
x = linspace(0,(params.N-1)*params.dz,params.N); % m
zoom_width = 10e-6; % m
frame_dec = 10;
nonlinear_f_plot(x, tspan, X', params, zoom_width, gen_video, plot_jacobian, frame_dec, 'nonlinear_f_test.avi');
