function X = generate_state_vector_nonlinear_wave_1D(E, dtE, P, dtP, p);
  % E_i = electric field at node i
  % dtE_i = time derivative of E_i
  % P_pi = pth-pole polarization response at node i
  % dtP_pi = time derivative of P_pi
  % X = state variables (ordering set by p.x_order)
  % params = parameters
  %  N: number of points in space (number of nodes)
  %  Lorentz: matrix of pole parameters
  %   (delta_chi_1, delta_1, omega_1),
  %   (delta_chi_2, delta_2, omega_2),
  %   ...
  %  x_order: 0 for default ordering above, 1 for ordering below
  %   0:  [E_1, dtE_1, P_11, dtP_11, P_21, dtP_21, P_31, dtP_31, ... E_2, dtE_2, P_12, dtP_12 ...]
  %   1:  [E_1, E_2, ... dtE_1, dtE_2, ... P_11, P_12, ... P_21, P_22, ... dtP_11, dtP_12, ... dtP_21, dtP_22]
  
  num_poles = size(params.Lorentz, 1);
  idx_step = 2*(1 + num_poles);
  X = zeros(params.N*idx_step, 1);
  if params.x_order == 0
    X(1:idx_step:end) = E;
    X(2:idx_step:end) = dtE;
    for p = 1:num_poles
      X(3+2*(p-1):idx_step:end) = P(p);
      X(4+2*(p-1):idx_step:end) = dtP(p);
    end
  else
    X(1:p.N) = E;
    X(p.N+1:2*p.N) = dtE;
    for p = 1:num_poles
      X(2*p.N+1 + p.N*(p-1):2*p.N + p.N*p) = P(p);
      X(p.N*(2+num_poles)+1 + p.N*(p-1):p.N*(2+num_poles) + p.N*p) = dtP(p);
    end
  end
  return X;
end
