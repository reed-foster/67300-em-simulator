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

  % material parameters
  params.eps_0 =  8.85e-12; % F/m
  params.mu_0 = 1.26e-6; % N/A^2
  params.chi_2 = 41.7e-12; % m/V Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
  params.chi_3 = 1.5e-20; % m^2/V^2 https://onlinelibrary.wiley.com/doi/pdf/10.1002/pssb.202200453 
  params.Lorentz = [4 2*pi*1e10 2*pi*1e11; ...
                    2 2*pi*9e13 2*pi*3e14]; % arbitrary

  % simulation params
  params.x_order = 0;
  params.N = 200; % discretization
  params.dz = 100e-9; % m
  params.dt = params.dz/3e8; % s
  params.tf = 300e-15; % s

  % source parameters
  params.omega_J = 2*pi*3e8/1.55e-6; % angular frequency for 1.55um
  params.t0_J = 2*2*pi/params.omega_J; % s
  params.ampl_J = 5e7/params.dz*params.omega_J; % arbitrary
end
