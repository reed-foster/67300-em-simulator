function [A, B] = Linearize_f(eval_f, x_0, p, eval_u, u)

N = length(x_0);

J_f = JacobianCalculation(eval_f, x_0, eps, N);
J_u = JacobianCalculation(eval_u, x_0, eps, N);

K_0 = eval_f(x_0) - J_f*x_0 - J_u*u;
B = [K_0 J_u];
A = J_f;

end