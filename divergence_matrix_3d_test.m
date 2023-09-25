dbg_plot = true; % plot 2D contourf of divergence of vector field

% test for divergence_matrix_2d method
[x,y] = meshgrid(-100:2:100,-20:2:20);
Fx = 200 - (2*sin(pi*x.^2/(50)) + 2*y.^2);
Fy = 200 - (10*sin(pi*x.^2/(50)) + 0.1*y.^2);

D = divergence(x,y,Fx,Fy);

% divide by twice the step size
div = reshape(divergence_matrix_2d(101,21)*[reshape(Fx',[2121 1]); reshape(Fy',[2121 1])], [101 21])' / 4;

trim = 2; % trim boundaries, since edge calculation is slightly different from matlab's implementation

D_trim = D(trim:end-trim,trim:end-trim);
div_trim = div(trim:end-trim,trim:end-trim);

err = norm(D_trim - div_trim);

disp(sprintf("error = %d", err));

if (dbg_plot)
  close all;
  subplot(2,1,1);
  contourf(x(trim:end-trim,trim:end-trim),y(trim:end-trim,trim:end-trim),D(trim:end-trim,trim:end-trim));
  
  subplot(2,1,2);
  contourf(x(trim:end-trim,trim:end-trim),y(trim:end-trim,trim:end-trim),div(trim:end-trim,trim:end-trim));
end 
