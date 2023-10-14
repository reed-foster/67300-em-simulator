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
  P = zeros(num_poles, p.N); % PxN
  dtP = zeros(num_poles, p.N); % PxN
  for p = 1:num_poles
    P(p) = X(3+2*(p-1):idx_step:end);
    dtP(p) = X(4+2*(p-1):idx_step:end);
  end

  % e.g. with 3 poles
  % X = [E_1, dtE_1, P_11, dtP_11, P_21, dtP_21, P_31, dtP_31, E_2, dtE_2, P_12, dtP_12 ...]
  %      1    2      3     4       5     6       7     8       9    10     11    12
  % idx_step = 2*(1 + 3) = 8;

  % eqn 1: rank reduction (d/dt(E) = dtE);
  dxdt(1:idx_step:end) = dtE;
  dtdtP = zeros(num_poles, p.N);
  for p = 1:num_poles
    % eqn 3: rank reduction (d/dt(P) = dtP);
    dxdt(3+2*(p-1):idx_step:end) = dtP(p);
    % eqn 4: dispersive polarization response (d/dt(dtP) = omega^2*eps0*delta_chi*E - omega^2*P - delta*dtP)
    dtdtP(p) = (p.P(p,3)^2*p.eps_0*p.P(p,1)).*E - (p.P(p,3)^2).*P(p) - (p.P(p,2)).*dtP(p);
    dxdt(4+2*(p-1):idx_step:end) = dtdtP(p);
  end
  % eqn 2: wave equation (d/dt(dtE) = 1/mu_0*(d/dz)(d/dz)E ...);
  dxdt(2:idx_step:end) = (1/p.mu_0*del2(E,dz^2) ...
                          - (2*p.eps_0).*(p.chi_2 + E.*(3*p.chi_3)).*(dtE.^2) ...
                          - sum(dtdtP, 1) ...
                          + dtJ) ./ (1 + (2*p.chi_2).*E + (3*p.chi_3).*(E.^2));
  return dxdt;
end
