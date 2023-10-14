function X = generate_states_nonlinear_wave_1D(E, dtE, P, dtP, p);
  % E_i = electric field at node i
  % dtE_i = time derivative of E_i
  % P_pi = pth-pole polarization response at node i
  % dtP_pi = time derivative of P_pi
  % X = [E_1, dtE_1, P_11, dtP_11, P_21, dtP_21, P_31, dtP_31, ... E_2, dtE_2, P_12, dtP_12 ...]
  % p: parameters
  %  N: number of points in space (number of nodes)
  %  P: matrix of pole parameters
  %   (delta_chi_1, delta_1, omega_1),
  %   (delta_chi_2, delta_2, omega_2),
  %   ...
  %  ...

  num_poles = size(p.P, 1);
  idx_step = 2*(1 + num_poles);
  X = zeros(p.N*idx_step, 1);
  X(1:idx_step:end) = E;
  X(2:idx_step:end) = dtE;
  for p = 1:num_poles
    X(3+2*(p-1):idx_step:end) = P(p);
    X(4+2*(p-1):idx_step:end) = dtP(p);
  end
  return X;
end
