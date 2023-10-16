function params = nonlinear_params;
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

  params.N = 4000; % discretization
  params.eps_0 =  8.85e-12; % F/m
  params.mu_0 = 1.26e-6; % N/A^2
  %params.chi_2 = 41.7e-12; % m/V Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
  params.chi_2 = 0;
  %params.chi_3 = 1.5e-20; % m^2/V^2 https://onlinelibrary.wiley.com/doi/pdf/10.1002/pssb.202200453 
  params.chi_3 = 0;
  %params.Lorentz = [4 2*pi*1e12 2*pi*1e13];
  params.Lorentz = [4 2*pi*9e18 2*pi*1e15];
  %params.Lorentz = [4 2*pi*1e10 2*pi*1e11; ...
  %                  2 2*pi*1e15 2*pi*1e16]; % arbitrary
  params.dz = 50e-9; % m
  params.x_order = 0;
end
