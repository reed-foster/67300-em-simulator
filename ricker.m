function ricker_wavelet = ricker(ampl, omega_p, N, dt, t0);
  tspan = linspace(0, (N-1)*dt, N);
  ricker_wavelet = ampl*(1 - 1/2*((tspan - t0)*omega_p).^2).*exp(-(0.5*(tspan-t0)*omega_p).^2);
end
