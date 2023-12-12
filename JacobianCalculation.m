function J_f = JacobianCalculation(f, x_0, eps, N, M)
  J_f = sparse(N, M);
  for i = 1:M
    e_vet = zeros(M, 1);
    e_vet(i) = 1;
    J_f(:, i) = 1./eps.*(f(x_0+e_vet.*eps)-f(x_0));
  end
end
