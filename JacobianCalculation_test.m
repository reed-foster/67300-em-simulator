clear
close all
clc

% parameters 
p.N = 25; % discretization
p.eps_0 =  8.85e-12; % F/m
p.mu_0 = 1.26e-6; % N/A^2
p.delta_x = 100e-9; % m

% nodal variables
E = zeros(p.N,1); % electric field
H = zeros(p.N,1); % magnetic field 

% source 
J = zeros(p.N,1);

% simulation setup and initial conditions
X0 = [sqrt(p.eps_0/p.mu_0)*gaussian_start(1,100,400,100,p.N)'; -gaussian_start(1,100,400,100,p.N)']; % V/m and A/m

% change logarithmicly the step
eps = logspace(10, -20, 100);
rel_err = zeros(length(eps), 1);

% analytical jacobian
TL = zeros(p.N, p.N);
TR = diag(ones(p.N-1, 1).*(1/(p.mu_0*2*p.delta_x)), 1) - diag(ones(p.N-1, 1).*(1/(p.mu_0*2*p.delta_x)), -1);
BL = diag(ones(p.N-1, 1).*(1/(p.eps_0*2*p.delta_x)), 1) - diag(ones(p.N-1, 1).*(1/(p.eps_0*2*p.delta_x)), -1);
BR = zeros(p.N, p.N);

JacobAnalytic = [TL, TR; BL, BR];

% compute Jacobian
for i = 1:length(eps)
    Jcalc = JacobianCalculation(@(X) evalf_freespace(X,J,p), X0, eps(i), 2*p.N);
    rel_err(i) = norm(Jcalc - JacobAnalytic)/norm(JacobAnalytic);
end

loglog(eps, rel_err, 'LineWidth', 1.5)
title('Relative error for Jacobian over $\epsilon$', 'Interpreter', 'latex')
xlabel('$\epsilon$', 'Interpreter', 'latex')
ylabel('$\frac{\|Jc - Ja\|}{\|Ja\|}$', 'Interpreter', 'latex')
grid on

%% print the Jacobian
custom_spy(JacobAnalytic)
print('spy_plot.png', '-dpng', '-r600');



