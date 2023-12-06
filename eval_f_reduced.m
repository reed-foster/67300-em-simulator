function dxdt = eval_f_reduced(x, p, u, A, B)
    dxdt = A*x + B*[1; u];
end
