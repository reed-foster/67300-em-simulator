function ricker_wavelet = ricker(t, ampl, omega_p, t0);
  ricker_wavelet = ampl*(1 - 1/2*((t - t0)*omega_p).^2).*exp(-(0.5*(t - t0)*omega_p).^2);
end
