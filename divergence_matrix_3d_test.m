% test for divergence_matrix_3d method
[x,y] = meshgrid(-8:2:8,-8:2:8);
Fx = 200 - (x.^2 + y.^2);
Fy = 200 - (x.^2 + y.^2);

D = divergence(x,y,Fx,Fy);

div = divergence_matrix_2d(9,9);

reshape(div,[9,9]);

subplots
contourf


spy(div)
