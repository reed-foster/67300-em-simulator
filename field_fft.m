function field_fft(E, z, t, p)
    color_order = lines(2);
    x = E;
    subplot(2, 1, 1);
    plot(z, x, 'linewidth', 2.5, 'Color', color_order(1, :));
    xlabel('Position (m)');
    ylabel('Amplitude');
    title("Signal in Time Domain at t = " + num2str(t) + " s")
    xlim([min(z), max(z)])
    ax = gca;
    ax.FontSize = 14;

    % Fourier transform
    subplot(2, 1, 2);
    Ts = p.dz;
    y = fft(x);
    fs = 1 / Ts;    
    n = length(x);
    fshift = (-n/2:n/2-1) * (fs/n);
    yshift = fftshift(y);
    % semilogy(fshift, abs(yshift), 'linewidth', 2, 'Color', color_order(2, :));
    plot(fshift, abs(yshift), 'linewidth', 2.5, 'Color', color_order(2, :));
    title("Signal in Wavevector Domain")
    xlabel('Wavevector (1/m)');
    ylabel('Amplitude');
    xlim([-max(fshift)/3, max(fshift)/3]);
    ax = gca;
    ax.FontSize = 14;
    % ylim([0, max(max(abs(fft(E))))])
    %ylim([min(min(abs(fft(E)))), max(max(abs(fft(E))))])
    % ylim([1e-6, 1e2]);
    
    % % Aggiungi una pausa tra i frame
    % pause(0.001);
    % 
    % frame = getframe(gcf);
    % M(k) = frame;
    % k = k+1;
end
