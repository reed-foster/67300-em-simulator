clear
close all

%% LOW FIELD
[p, newton_opts, trap_opts] = demo_initialize(5e-6, 501);

p.ampl_J = 1e4/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-16;                   % minimum timestep
p.tf = 42e-15;                  % final simulation time
p.source_type = "ricker";       % source type

% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = 2e-16;

% compute preconditioner
Jf_rick_low = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% trapezoidal method
[t_rick_low,X_rick_low] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf_rick_low, trap_opts, newton_opts);

%% HIGH FIELD
[p, newton_opts, trap_opts] = demo_initialize(5e-6, 501);

p.ampl_J = 3e7/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-16;                   % minimum timestep
p.tf = 42e-15;                  % final simulation time
p.source_type = "ricker";       % source type

% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = 2e-16;

% compute preconditioner
Jf_rick_high = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% trapezoidal method
[t_rick_high, X_rick_high] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf_rick_high, trap_opts, newton_opts);

%% LOW FIELD SINUSOIDAL

[p, newton_opts, trap_opts] = demo_initialize(20e-6, 501);

p.ampl_J = 1e4/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-17;                   % minimum timestep
p.tf = 60e-15;                  % final simulation time
p.source_type = "sinusoid";       % source type

% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = Inf;

tic
% compute preconditioner
Jf_sin_low = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% trapezoidal method
[t_sin_low, X_sin_low] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf_sin_low, trap_opts, newton_opts);
toc

%% FFT SINUSOIDAL
movie_fft(X_sin_low, t_sin_low, p, false)

%% HIGH FIELD SINUSOIDAL

[p, newton_opts, trap_opts] = demo_initialize(20e-6, 501);

p.ampl_J = 3e7/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-16;                   % minimum timestep
p.tf = 150e-15;                  % final simulation time
p.source_type = "sinusoid";       % source type

% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = Inf;

tic
% compute preconditioner
Jf_sin_high = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% trapezoidal method
[t_sin_high, X_sin_high] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf_sin_high, trap_opts, newton_opts);
toc

%% FFT SINUSOIDAL
[p, ~, ~] = demo_initialize(20e-6, 501);
movie_fft(X_sin_high, t_sin_high, p, false)

%% COLORMAP
[p, ~, ~] = demo_initialize(20e-6, 501);
space_colormap(X_sin_high, t_sin_high, p)

%% 5000 NODES
[p, newton_opts, trap_opts] = demo_initialize(2e-5, 5001);

p.ampl_J = 3e7/p.dz*p.omega_J;  % source amplitude
p.dt = 2e-17;                   % minimum timestep
p.tf = 30e-15;                  % final simulation time
p.source_type = "ricker";       % source type

% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = 2e-16;

% compute preconditioner
Jf_long = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% trapezoidal method
[t_long,X_long] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf_long, trap_opts, newton_opts);

