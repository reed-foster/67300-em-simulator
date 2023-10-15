function gauss_wavelet = gaussian_start(ampl, width, pos, sp_wave, N, dz)
    % ampl = amplitude of the gaussian [a.u.]
    % width = standard deviation of the gaussian [m]
    % pos = mean of the gaussian [m]
    % sp_wave = wavelength of the sinusoidal modulation [m]
    % N = number of points
    % dz = distance between points [m]
    f = @(x, mu, std) ampl.*exp(-(x-mu).^2./(2*std.^2)).*cos(2*pi*(x-pos)/sp_wave);
    x = linspace(0,(N-1)*dz,N);
    gauss_wavelet = f(x, pos, width);
end
