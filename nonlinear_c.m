function [E, D, dtE, dtD, P, dtP] = nonlinear_c(X,params);
  % get the output from the state X
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
  
  [E, dtE, P, dtP] = nonlinear_split_X(X,params);
  D = params.eps_0*(E + params.chi_2*E.^2 + params.chi_3*E.^3) + sum(P,1)';
  dtD = params.eps_0*(dtE + 2*params.chi_2*E.*dtE + 3*params.chi_3*E.^2.*dtE) + sum(dtP,1)';
end
