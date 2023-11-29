clear all
close all

p = simulation_params();

% simulation setup and initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);

p.ampl_J = 1e6/p.dz*p.omega_J;
p.dt = 1e-18;

newton_opts.err_f = Inf;
newton_opts.err_dx = Inf;
newton_opts.err_rel = 1e-2;
newton_opts.max_iter = 20;
newton_opts.matrix_free = true;
newton_opts.err_gcr = 1e-2; % relative error for GCR residual inside Newton
newton_opts.eps_fd = 1e-7; % relative perturbation for Jacobian
newton_opts.save_intermediate = false;

trap_opts.save_intermediate = true;
trap_opts.visualize = 1000;

tic;
[t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, trap_opts, newton_opts);
toc;
