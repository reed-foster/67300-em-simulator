function custom_spy(M)
    Pos = M > 0;
    Neg = M < 0;
    
    % Create a custom spy plot
    figure;
    spy(Pos, 'r')%'color', Red); % Use red color for negative values
    hold on;
    spy(Neg, 'b')%'color', Blu); % Use blue color for positive values
    title('Matrix sparsity');
        
    % Modify the appearance of the plot if necessary
    axis square;
end
