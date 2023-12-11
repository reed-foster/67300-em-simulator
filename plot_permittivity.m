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
  end
  yyaxis left;
  semilogx(omega/(2*pi), real(eps_r), LineWidth=2);
  ylabel("real relative permittivity \epsilon'");
  yyaxis right;
  semilogx(omega/(2*pi), imag(eps_r), LineWidth=2);
  ylabel("imaginary relative permittivity \epsilon''");
  xline([1 3 5 7 9]*params.omega_J/(2*pi), '--');
  legend("\epsilon'", "\epsilon''", "\omega_J");
  xlabel("freq [Hz]");
  xlim([min(omega/(2*pi)) max(omega/(2*pi))]);
  title("Lorentz 2-oscillator complex permittivity");
  ax = gca;
  ax.FontSize = 18;
end
