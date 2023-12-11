function p = simulation_params;
  % p = parameters
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
  p.eps_0 =  8.85e-12; % F/m
  p.mu_0 = 1.26e-6; % N/A^2
  p.chi_2 = 41.7e-12; % m/V Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
  p.chi_3 = 1.5e-20; % m^2/V^2 https://onlinelibrary.wiley.com/doi/pdf/10.1002/pssb.202200453 
  p.Lorentz = [4 2*pi*1e10 2*pi*1e11; ...
               2 2*pi*3e14 2*pi*3e15]; % arbitrary

  % simulation parameters
  p.x_order = 0;
  p.N = 3001; % discretization
  p.dz = 10e-9; % m
  p.dt = p.dz/3e8; % s
  p.tf = 20e-15; % s

  % source parameters
  p.omega_J = 2*pi*3e8/1.55e-6; % angular frequency for 1.55um
  p.t0_J = 2*2*pi/p.omega_J; % s, sets delay for ricker wavelet, also sets turnon time of source
  p.ampl_J = 5e7/p.dz*p.omega_J; % arbitrary
  %p.ampl_J = 5e-3/p.dz*p.omega_J; % arbitrary
  p.source_type = "ricker"; % ricker or sinusoid

  % precompute unit scaling
  p.E_scale = p.eps_0;
  p.P_scale = 1;
  p.dtE_scale = p.E_scale*p.dt*100;
  p.dtP_scale = p.P_scale*p.dt*100;
  e = ones(p.N,1);
  p.X_scale = generate_X(e*p.E_scale, e*p.dtE_scale, e*p.P_scale, e*p.dtP_scale, p);
end
