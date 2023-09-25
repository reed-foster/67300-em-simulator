dbg_plot = false; % plot 2D contourf of divergence of vector field

% test for divergence_matrix_2d method
[x,y] = meshgrid(-100:2:100,-20:2:20);
Fx = 200 - (2*sin(pi*x.^2/(50)) + 2*y.^2);
Fy = 200 - (10*sin(pi*x.^2/(50)) + 0.1*y.^2);

% matlab implementation of divergence (does a central difference in the middle and left/right-handed difference on the boundaries)
D = divergence(x,y,Fx,Fy);

% our finite difference implementation
% divide by twice the step size
raw = divergence_matrix_2d(101,21)*[reshape(Fx',[2121 1]); reshape(Fy',[2121 1])];
div = reshape(raw, [101 21])' / 4;

% trim boundaries, since edge calculation is slightly different from matlab's implementation
% we don't really care too much about exactly what the divergence calculation will be at the boundary, since the field should be close to zero when we use a PML
% if we add periodic boundary conditions, then we can change the matrix so that the difference terms which wrap around a boundary are included.
trim = 2;

D_trim = D(trim:end-trim,trim:end-trim);
div_trim = div(trim:end-trim,trim:end-trim);

err = norm(D_trim - div_trim);

disp(sprintf("2D error = %d", err));

if (dbg_plot)
  close all;
  subplot(2,1,1);
  contourf(x(trim:end-trim,trim:end-trim),y(trim:end-trim,trim:end-trim),D_trim);
  
  subplot(2,1,2);
  contourf(x(trim:end-trim,trim:end-trim),y(trim:end-trim,trim:end-trim),div_trim);
end 

% test for divergence_matrix_2d method
[x,y,z] = meshgrid(-10:2:10,-20:2:20,-70:2:70);
Fx = 200 - (2*sin(pi*x.^2/(50)) + 2*y.^2 + 10*exp(-z.^2/100));
Fy = 200 - (10*sin(pi*x.^2/(50)) + 0.1*y.^2 + 2*exp(-z.^2/10));
Fz = 200 - (0.1*sin(pi*x.^2/(50)) + 4*y.^2 + 0.1*exp(-z.^2/800));

% matlab implementation of divergence (does a central difference in the middle and left/right-handed difference on the boundaries)
D = divergence(x,y,z,Fx,Fy,Fz);

% our finite difference implementation
% divide by twice the step size
raw = divergence_matrix_3d(11,21,71)*[reshape(permute(Fx,[3 2 1]),[16401 1]); reshape(permute(Fy,[3 2 1]),[16401 1]); reshape(permute(Fz,[3 2 1]),[16401 1])];
div = permute(reshape(raw, [11 21 71]), [3 2 1]) / 4;

% trim boundaries, since edge calculation is slightly different from matlab's implementation
% we don't really care too much about exactly what the divergence calculation will be at the boundary, since the field should be close to zero when we use a PML
% if we add periodic boundary conditions, then we can change the matrix so that the difference terms which wrap around a boundary are included.
trim = 2;

D_trim = D(trim:end-trim,trim:end-trim,trim:end-trim);
div_trim = div(trim:end-trim,trim:end-trim,trim:end-trim);

err = norm(D_trim - div_trim);

disp(sprintf("3D error = %d", err));
