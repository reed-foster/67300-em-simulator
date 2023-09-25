function div = divergence_matrix_3d(L,M,N)
  % returns the 3D central finite-difference approximation of the divergence of a vector field.
  % the matrix transforms a 1D state vector ordered first by increasing spatial coordinate along the x,y,z axes (respectively),
  % then by vector direction e.g. V_x(0,0,0) V_x(1,0,0) ... V_x(0,1,0) ... V_y(0,0,0) ...
  % L is the number of points in space in the x direction, M the number of points in the y direction, and N in the z direction
  % use speye instead of eye to reduce matrix size (in memory) and computation time
  div = [kron(speye(N*M),divergence_matrix_1d(L)) kron(kron(speye(N),divergence_matrix_1d(M)),speye(L)) kron(divergence_matrix_1d(N),speye(M*L))];
end
