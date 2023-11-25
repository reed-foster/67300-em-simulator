function [E, dtE, P, dtP] = split_X(X, params);
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
  dxdt = zeros(params.N*idx_step, 1);
  % get the offset into the state variable
  if params.x_order == 0
    E = X(1:idx_step:end); % Nx1
    dtE = X(2:idx_step:end); % Nx1
    P = zeros(num_poles, params.N); % PxN
    dtP = zeros(num_poles, params.N); % PxN
    for p = 1:num_poles
      P(p,:) = X(3+2*(p-1):idx_step:end);
      dtP(p,:) = X(4+2*(p-1):idx_step:end);
    end
  else
    E = X(1:params.N);
    dtE = X(params.N+1:2*params.N);
    for p = 1:num_poles
      P(p,:) = X(2*params.N+1 + params.N*(p-1):2*params.N + params.N*p);
      dtP(p,:) = X(params.N*(2+num_poles)+1 + params.N*(p-1):params.N*(2+num_poles) + params.N*p);
    end
  end
end
