clear all
close all

p = params();

% source 
p.tsteps = round(p.tf/p.dt);
p.tvec = linspace(0,(p.tsteps-1)*p.dt,p.tsteps); % s

% simulation setup and initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);

options = odeset('RelTol',1e-3,'AbsTol',1e-12);
eval_f = @(t,X) eval_f(X,eval_u(t,p),p);
[t,X] = ode45(eval_f, p.tvec, X0, options);

gen_video = true;
plot_jacobian = false;
xvec = linspace(0,(p.N-1)*p.dz,p.N); % m
zoom_width = 10e-6; % m
frame_dec = 12;
visualize_state(xvec, p.tvec, X', p, zoom_width, gen_video, plot_jacobian, frame_dec, 'test.avi');
