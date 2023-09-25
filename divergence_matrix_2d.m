function div = divergence_matrix_2d(L,M)
  % returns the 2D central finite-difference approximation of the divergence of a vector field.
  % the matrix transforms a 1D state vector ordered first by increasing spatial coordinate along the x and y axes (respectively),
  % then by vector direction e.g. V_x(0,0) V_x(1,0) ... V_x(0,1) ... V_y(0,0) ...
  % L is the number of points in space in the x direction, M the number of points in the y direction
  div = [kron(eye(M),divergence_matrix_1d(L)) kron(divergence_matrix_1d(M), eye(L))];
end
