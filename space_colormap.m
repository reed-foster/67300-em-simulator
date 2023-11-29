function space_colormap(X, t, p)

    % Change when vectorize split_X
    for i = 1:size(X, 2)
        [E(:, i), ~, ~, ~] = split_X(X(:, i), p);
    end
    
    % Change it save z somewhere
    z_temp = 0:p.dz:(p.dz*size(E, 1)-p.dz);
    z = z_temp - mean(z_temp);
    
    [T, Z] = meshgrid(t, z);
    figure;
    surf(T, Z, E, 'EdgeColor','None');
    view(2);
    xlabel('Time (s)');
    ylabel('Position (m)');
    zlabel('Valore E');
    title('Electric field space and time evolution');
    colormap(brewermap([],'RdBu'))
    axis vis3d
end

