function dxdt = evalf_nonlinear_wave_1D(X,dtJ,params)
  % E_i = electric field at node i
  % dtE_i = time derivative of E_i
  % P_pi = pth-pole polarization response at node i
  % dtP_pi = time derivative of P_pi
  % X = state variables (ordering set by p.x_order)
  % dtJ = input (time derivative of source term)
  % params = parameters
  %  N: number of points in space (number of nodes)
  %  dz: distance between points in space
  %  Lorentz: matrix of pole parameters
  %   (delta_chi_1, delta_1, omega_1),
  %   (delta_chi_2, delta_2, omega_2),
  %   ...
  %  mu_0: vaccuum permeability
  %  eps_0: vaccuum permittivity
  %  chi_2: second order (instantaneous) nonlinearity
  %  chi_3: third order (instantaneous) nonlinearity
  %  x_order: 0 for default ordering above, 1 for ordering below
  %   0:  [E_1, dtE_1, P_11, dtP_11, P_21, dtP_21, P_31, dtP_31, ... E_2, dtE_2, P_12, dtP_12 ...]
  %   1:  [E_1, E_2, ... dtE_1, dtE_2, ... P_11, P_12, ... P_21, P_22, ... dtP_11, dtP_12, ... dtP_21, dtP_22]

  num_poles = size(params.Lorentz, 1);
  idx_step = 2*(1 + num_poles);
  dxdt = zeros(params.N*idx_step, 1);
  % split up state variable X so we can do vectorized computation of dxdt
  [E, dtE, P, dtP] = split_state_vector_nonlinear_wave_1D(X, params);

  % eqn 1: rank reduction (d/dt(E) = dtE);
  dxdt(1:idx_step:end) = dtE;
  dtdtP = zeros(num_poles, params.N);
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
