function space_colormap(X, t, p)
    n = 1;
    
    % Change when vectorize split_X
    for i = 1:size(X, 2)
        [E(:, i), ~, ~, ~] = split_X(X(:, i)./p.X_scale, p);
    end
    
    % Change it save z somewhere
    z_temp = 0:p.dz:(p.dz*size(E, 1)-p.dz);
    z = z_temp - mean(z_temp);
    
    [T, Z] = meshgrid(t(1:n:end)*1e15, z(1:n:end)*1e6);
    figure;
    surf(T, Z, E(1:n:end,1:n:end), 'EdgeColor','None');
    view(2);
    xlabel('Time [fs]');
    ylabel('Position [um]');
    zlabel('Valore E');
    h = colorbar;
    h.Label.String = 'E-field [V/m]';
    m_E = max(abs(E),[],'all')*1.1;
    caxis([-m_E m_E]);
    title('Electric field space and time evolution');
    colormap(brewermap([],'RdBu'))
    ax = gca;
    ax.FontSize = 18;
end

