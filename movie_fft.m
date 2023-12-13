function movie_fft(X, t, p, movie)

    % Change when vectorize split_X
    for i = 1:size(X, 2)
        [E(:, i), ~, ~, ~] = split_X(X(:, i), p);
    end
    
    % Change it save z somewhere
    z_temp = 0:p.dz:(p.dz*size(E, 1)-p.dz);
    z = z_temp - mean(z_temp);
    
    figure;
    k = 1;
    for i = 1:100:length(E(1, :))
        
        fft_fig = figure(1);
        set(fft_fig, 'resize', 'off', 'Position', [100 100 1380 820]);
        field_fft(E(:, i), z, t(i), p)
        pause(0.1)
    
        frame = getframe(gcf);
        M(k) = frame;
        k = k+1;
    end

    if movie == true
        Prova = M(floor(1:0.3:end));
    
        % Create a VideoWriter object
        writerObj = VideoWriter('TestMovie.mp4', 'MPEG-4');
        open(writerObj);
        
        % Write frames to the video file
        for i = 1:length(Prova)
            writeVideo(writerObj, Prova(i).cdata);
        end
        
        % Close the video file
        close(writerObj);
    end
end