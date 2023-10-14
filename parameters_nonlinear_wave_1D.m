function p = parameters_nonlinear_wave_1D()

  %parameters 
  p.N = 1500; % discretization
  p.eps_0 =  8.85e-12; % F/m
  p.mu_0 = 1.26e-6; % N/A^2
  p.chi_2 = 41.7e-12; % m/V Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
  p.chi_3 = 1.5e-20; % m^2/V^2 https://onlinelibrary.wiley.com/doi/pdf/10.1002/pssb.202200453 
  p.P = [4 2*pi*1.2e12 2*pi*2.9e14]; % arbitrary
  p.dz = 100e-9; % m
  
  return p
end
