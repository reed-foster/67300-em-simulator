clear
close all
clc

A = 1;
B = 10;
C = 3;
D = 2;
E = 4;
F = 5;

f = @(x) [A*x(1) + B*x(2) + x(3); C*x(1) + D*x(2) + x(3); E*x(1) + F*x(2) + x(3)];

eps = 0.001;

Start = [f([1, 0, 0]), f([0, 1, 0]), f([0, 0, 1])]
Jacob = JacobianCalculation(@(x) f(x), [1, 1, 1], eps, 3, 3)
