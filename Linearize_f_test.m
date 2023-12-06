% Code to check Linearize_f

x_0 = [1; 2];
u_0 = [3; 4];
p.N = length(x_0);

x = [1; 2];
u = [3; 4];

f = @(x, p, u) [3*x(1)^2; 6 + 2*x(2)^2 + u(1)^3];
Jfx = [6*x_0(1), 0; 0, 4*x_0(2)];
Jfu = [0, 0; 3*u_0(1)^2, 0];
K0 = f(x_0, p, u_0) - Jfx*x_0 - Jfu*u_0;
Aex = Jfx;
Bex = [K0, Jfu];

eps = 1e-6;

[A, B] = Linearize_f(f, x_0, p, u_0, eps);

state = A*x + B*[1; u];

A
Aex

state
f(x, p, u)
eval_f_reduced(x, p, u, A, B)


