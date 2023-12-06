clear all
close all
clc

f = @(x) x^3;

x_0 = 1;
u_0 = 2;

eval_f = f;
eval_u = @(x) u_0;

eps = 1e-6;
N = length(x_0);

[A, B] = Linearize_f(eval_f, x_0, eps, eval_u, u_0);

df_dx_numeric = (f(x_0 + eps) - f(x_0 - eps)) / (2 * eps);

disp(A);
disp(df_dx_numeric);
