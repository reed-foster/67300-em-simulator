function J_f = JacobianCalculation(f, x_0, eps, N)
    J_f = zeros(N, N);
    for i = 1:N
        e_vet = zeros(N, 1);
        e_vet(i) = 1;
        J_f(:, i) = 1./eps.*(f(x_0+e_vet.*eps)-f(x_0));
    end
end