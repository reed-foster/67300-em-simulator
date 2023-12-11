function T_inv = trap_newton_gcr_preconditioner(Jf0, p);
  % computes a left preconditioner matrix inv(T) for use with GCR
  % usage
  % trap_newton_gcr_preconditioner(Jf0, p)
  % Jf0: Jacobian of system under zero field
  % p: parameters of system

  Jf0_trap = speye(size(Jf0, 1)) - Jf0*p.dt/2;
  N_nodal_vars = 2+2*size(p.Lorentz,1);
  % pick a node away from the edge so we don't have the Dirichlet boundary condition
  T_node = Jf0_trap(N_nodal_vars+1:2*N_nodal_vars,N_nodal_vars+1:2*N_nodal_vars);
  T_inv = kron(speye(p.N), sparse(T_node)\speye(N_nodal_vars));
end
