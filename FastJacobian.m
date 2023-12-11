function Jf = FastJacobian(f, x0, p, dx, N);
  N_sub = 2 + 2*size(p.Lorentz,1); % E, dtE, P1, dtP1, P2, dtP2, etc.
  Jf = sparse(N_sub*N, N_sub*N);
  e_vet = sparse(N*N_sub, 1);
  first = sparse(N_sub,N_sub);
  middle = sparse(N_sub,N_sub);
  last = sparse(N_sub,N_sub);
  lower = sparse(N_sub,N_sub);
  upper = sparse(N_sub,N_sub);
  for i = [1:2*N_sub (N-1)*N_sub+1:N*N_sub]
    if i > 1
      e_vet(i-1) = 0;
    end
    e_vet(i) = 1;
    dxdt = 1./dx.*(f(x0 + e_vet.*dx)-f(x0));
    if (i <= N_sub)
      first(:, i) = dxdt(1:N_sub);
      lower(:, i) = dxdt(N_sub+1:2*N_sub);
    elseif (i <= 2*N_sub)
      upper(:, i-N_sub) = dxdt(1:N_sub);
      middle(:, i-N_sub) = dxdt(N_sub+1:2*N_sub);
    else
      last(:, i-(N-1)*N_sub) = dxdt(end-N_sub+1:end);
    end
  end
  Jf(1:N_sub,1:N_sub) = first;
  Jf(N_sub+1:end-N_sub,N_sub+1:end-N_sub) = kron(speye(N - 2), middle);
  Jf(end-N_sub+1:end,end-N_sub+1:end) = last;
  Jf = Jf + kron(spdiags(ones(N,1), -1, N, N), lower) + kron(spdiags(ones(N,1), 1, N, N), upper);
end
