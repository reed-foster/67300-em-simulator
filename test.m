clear all
close all

p = simulation_params();

% simulation setup and initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);

p.ampl_J = 1e2/p.dz*p.omega_J;
p.dt = 1e-16;
p.tf = 20e-15;

newton_opts.err_f = Inf;
newton_opts.err_dx = Inf;
newton_opts.err_rel = 1e-2;
newton_opts.max_iter = 20;
newton_opts.matrix_free = true;
newton_opts.err_gcr = 1e-2; % relative error for GCR residual inside Newton
newton_opts.eps_fd = 1e-7; % relative perturbation for Jacobian
newton_opts.preconditioner = true;

trap_opts.save_intermediate = true;
trap_opts.visualize_dt = 2e-16;
%trap_opts.visualize_dt = Inf;
trap_opts.adaptive_timestep = false;
trap_opts.linear_only = false;
trap_opts.print_debug = false;

[Jf0, ~] = Linearize_f(@eval_f, X0, p, zeros(p.N,1), 1e-6);

tic;
[t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf0, trap_opts, newton_opts);
toc;

space_colormap(X, t, p);
figure;
plot_permittivity(p);
