clear
close all

params = nonlinear_params();
params.N = 1000;

% source 
dt = params.dz/3e8; %50e-18; % s
tf = 500e-15; % s
tsteps = round(tf/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s
J = zeros(params.N,tsteps);
dtJ = zeros(params.N,tsteps);

omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
dtJ_ricker = 2e8/params.dz*omega_J*(1-((tspan-2*pi/omega_J)*omega_J).^2).*exp(-(tspan-2*pi/omega_J).^2/(2*(1/omega_J)^2));
dtJ(round(params.N/2),:) = dtJ_ricker;

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.dz,params.N); % m

% E0 = zeros(params.N,1);
% dtE0 = zeros(params.N, 1);
% P0 = zeros(size(params.Lorentz,1), params.N);
% dtP0 = zeros(size(params.Lorentz,1), params.N);

E0 = ones(params.N,1);
dtE0 = 2*ones(params.N, 1);
P0 = 3*ones(size(params.Lorentz,1), params.N);
dtP0 = 4*ones(size(params.Lorentz,1), params.N);

X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

J_f = JacobianCalculation(@(X) nonlinear_f(X,dtJ(:,1),params), X0, 1e-4, 6*params.N);
custom_spy(J_f)

print('spy_plot.png', '-dpng', '-r600');
