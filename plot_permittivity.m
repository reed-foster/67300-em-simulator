function plot_permittivity(p)
  % p contains the following fields:
  % p.chi_1
  % p.omega_chi_1
  % p.delta_chi_1
  omega = 2*pi*logspace(6, 18, 500); % Hz
  eps_r = 1 + (p.chi_1*p.omega_chi_1^2)./(p.omega_chi_1^2-omega.^2+omega.*(1j*p.delta_chi_1));
 
  yyaxis left;
  semilogx(omega/(2*pi), real(eps_r));
  yyaxis right;
  semilogx(omega/(2*pi), imag(eps_r));
  legend("\epsilon'", "\epsilon''");
  xlabel("freq [Hz]");
end
