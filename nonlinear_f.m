function dxdt = nonlinear_f(X,dtJ,params)
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
  [E, dtE, P, dtP] = nonlinear_split_X(X, params);

  % eqn 1: rank reduction (d/dt(E) = dtE);
  dtdtP = zeros(num_poles, params.N);
  dchi_p = params.Lorentz(:,1);
  delta_p = params.Lorentz(:,2);
  omega_p = params.Lorentz(:,3);
  for p = 1:num_poles
    % eqn 4: dispersive polarization response (d/dt(dtP) = omega^2*eps0*delta_chi*E - omega^2*P - delta*dtP)
    dtdtP(p,:) = omega_p(p)^2*params.eps_0*dchi_p(p)*E - omega_p(p)^2*P(p) - delta_p(p)*dtP(p);
  end
  %disp(max(sum(dtdtP, 1), [], 'all'));
  % eqn 2: wave equation (d/dt(dtE) = 1/mu_0*(d/dz)(d/dz)E ...);
  dtdtE = 1/params.mu_0*del2(E,params.dz);
  %dtdtE = dtdtE + dtJ;
  %dtdtE = dtdtE - (2*params.eps_0).*(params.chi_2 + E.*(3*params.chi_3)).*(dtE.^2);
  dtdtE = dtdtE - sum(dtdtP, 1)' + dtJ;
  dtdtE = dtdtE ./ (params.eps_0*(1 + (2*params.chi_2).*E + (3*params.chi_3).*(E.^2)));
  dxdt = nonlinear_generate_X(dtE, dtdtE, dtP, dtdtP, params);
end
