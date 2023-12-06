function [A, B] = Linearize_f(eval_f, x_0, p, u_0, eps)

J_f = JacobianCalculation_Rect(@(x) feval(eval_f, x, p, u_0), x_0, eps, length(eval_f(x_0, p, u_0)), length(x_0));
J_u = JacobianCalculation_Rect(@(u) feval(eval_f, x_0, p, u), u_0, eps, length(eval_f(x_0, p, u_0)), length(u_0));

K_0 = eval_f(x_0, p, u_0) - J_f*x_0 - J_u*u_0;
B = [K_0 J_u];
A = J_f;

A = sparse(A);
B = sparse(B);

end