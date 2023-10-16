function plot_permittivity(params)
  % p contains the following fields:
  % p.chi_1
  % p.omega_chi_1
  % p.delta_chi_1
  omega = 2*pi*logspace(6, 18, 500); % Hz
  dchi_p = params.Lorentz(:,1);
  delta_p = params.Lorentz(:,2);
  omega_p = params.Lorentz(:,3);
  eps_r = 1;
  for p = 1:size(params.Lorentz,1)
    eps_r = eps_r + (dchi_p(p)*omega_p(p)^2)./(omega_p(p)^2-omega.^2-omega.*(1j*delta_p(p)));
 
  yyaxis left;
  semilogx(omega/(2*pi), real(eps_r));
  yyaxis right;
  semilogx(omega/(2*pi), imag(eps_r));
  legend("\epsilon'", "\epsilon''");
  xlabel("freq [Hz]");
end
