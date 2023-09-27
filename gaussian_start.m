function gauss_wavelet = gaussian_start(ampl, width, pos, sp_wave, N)
    % ampl = amplitude of the gaussian
    % width = standard deviation of the gaussian
    % pos = mean of the gaussian
    % sp_wave = wavelength of the sinusoidal modulation
    % N = number of points
    f = @(x, mu, std) ampl.*exp(-(x-mu).^2./(2*std.^2)).*sin(2*pi*x/sp_wave);
    x = 1:N;
    gauss_wavelet = f(x, pos, width);
end