clear
close all

params = nonlinear_params();
params.x_order = 0;
params.N = 100;

% source 
dt = params.dz/3e8; %50e-18; % s
tf = 40e-15; % s
tsteps = round(tf/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s
J = zeros(params.N,tsteps);
dtJ = zeros(params.N,tsteps);

omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
t0_J = 2*2*pi/omega_J;
dtJ(round(params.N/2),:) = ricker(5e7/params.dz*omega_J, omega_J, tsteps, dt, t0_J);

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.dz,params.N); % m

E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz,1), params.N);
dtP0 = zeros(size(params.Lorentz,1), params.N);

X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

options = odeset('RelTol',1e-3,'AbsTol',1e-12);

eval_f = @(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params);
[t,X] = ode45(eval_f, tspan, X0, options);

% change logarithmicly the step
eps = logspace(10, -10, 100);
rel_err = zeros(length(eps), 2);

% try for zero (initial condition) and nonzero (final state after perturbing material with free current)
t_idx = [1 size(X,1)];
for j = 1:2;
    Jbest = JacobianCalculation(@(X) nonlinear_f(X,0,params), X(t_idx(j),:)', 1e3, size(X,2), size(X,2));
    % compute Jacobian
    for i = 1:length(eps)
        Jcalc = JacobianCalculation(@(X) nonlinear_f(X,dtJ(:,end),params), X(t_idx(j),:)', eps(i), size(X,2), size(X,2));
        rel_err(i,j) = norm(Jcalc - Jbest)/norm(Jbest);
    end
end

loglog(eps, rel_err(:,1), 'LineWidth', 1.5);
hold on;
loglog(eps, rel_err(:,2), 'LineWidth', 1.5);
legend("zero-bias", "non-zero pulse state");
title('Relative error for Jacobian vs $\epsilon$', 'Interpreter', 'latex')
xlabel('$\epsilon$', 'Interpreter', 'latex')
ylabel('$\frac{\|Jc - J_0\|}{\|J_0\|}$', 'Interpreter', 'latex')

custom_spy(JacobianCalculation(@(X) nonlinear_f(X,0,params), X0, 1e3, size(X0,1), size(X0,1)));
title("Sparsity with new ordering");

params.x_order = 1;
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);
custom_spy(JacobianCalculation(@(X) nonlinear_f(X,0,params), X0, 1e3, size(X0,1), size(X0,1)));
title("Sparsity with old ordering");
