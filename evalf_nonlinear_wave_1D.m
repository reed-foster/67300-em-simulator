function dxdt = evalf_nonlinear_wave_1D(X,dtJ,p)
  % E_i = electric field at node i
  % dtE_i = time derivative of E_i
  % P_pi = pth-pole polarization response at node i
  % dtP_pi = time derivative of P_pi
  % X = [E_1, dtE_1, P_11, dtP_11, P_21, dtP_21, P_31, dtP_31, ... E_2, dtE_2, P_12, dtP_12 ...]
  % dtJ: input (time derivative of source term)
  % p: parameters
  %  N: number of points in space (number of nodes)
  %  dz: distance between points in space
  %  P: matrix of pole parameters
  %   (delta_chi_1, delta_1, omega_1),
  %   (delta_chi_2, delta_2, omega_2),
  %   ...
  %  mu_0: vaccuum permeability
  %  eps_0: vaccuum permittivity
  %  chi_2: second order (instantaneous) nonlinearity
  %  chi_3: third order (instantaneous) nonlinearity

  num_poles = size(p.P, 1);
  idx_step = 2*(1 + num_poles);
  dxdt = zeros(p.N*idx_step, 1);
  % get the offset into the state variable
  E = X(1:idx_step:end); % Nx1
  dtE = X(2:idx_step:end); % Nx1
  P = zeros(p.N, num_poles); % NxP
  dtP = zeros(p.N, num_poles); % NxP
  for i = 1:p.N
    for p = 1:p.P
      P(i,p) = X(3+idx_step*(i-1)+2*(p-1));
      dtP(i,p) = X(3+idx_step*(i-1)+2*(p-1)+1);
    end
  end

  for i = 1:p.N
    % eqn 1: rank reduction
    dxdt(i*idx_step+E_idx_offset) = dtE;
    % eqn 2: wave equation
    dxdt(i*idx_step+dtE_idx_offset) = (1/p.mu_0*laplacian_1D(i,idx_step,E_idx_offset,X,p.N,p.dz) ...
                                        - 2*p.eps_0*(p.chi_2 + 3*p.chi_3*X(;
  end
end
