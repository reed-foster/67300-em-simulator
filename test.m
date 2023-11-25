clear all
close all

p = params();

% simulation setup and initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);

newton_opts.err_f = Inf;
newton_opts.err_dx = 1e-8;
newton_opts.err_rel = Inf;
newton_opts.max_iter = 20;
newton_opts.save_intermediate = false;

trap_opts.save_intermediate = true;
trap_opts.visualize = 2;

eps_J = 1e-4;
Jf = @(x,p,u) JacobianCalculation(@(x) eval_f(x,p,u), x, eps_J, size(x,1));
[t,X] = trapezoid(@eval_f, Jf, p, @(t) eval_u(t,p), X0, p.tf, p.dt, trap_opts, newton_opts);
