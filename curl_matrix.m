function C = curl_matrix(L,M,N)
  % computes the curl matrix using sparse matrices
  % returns the central finite difference approximation of the curl of a vector field in R3
  % the matrix transforms a 1D state vector ordered first by increasing spatial coordinate along the x,y,z axes (respectively),
  % then by vector direction e.g. V_x(0,0,0) V_x(1,0,0) ... V_x(0,1,0) ... V_y(0,0,0) ...
  % L is the number of points in space in the x direction, M the number of points in the y direction, and N in the z direction
  curl_matrix = [0 -1 1; 1 0 -1; -1 1 0];
  P = L*M*N;
  ddx = kron(speye(N*M),divergence_matrix_1d(L));
  ddy = kron(kron(speye(N),divergence_matrix_1d(M)),speye(L));
  ddz = kron(divergence_matrix_1d(N),speye(M*L));
  ddx_stamp = [0 0 0; 0 0 -1; 0 1 0];
  ddy_stamp = [0 0 1; 0 0 0; -1 0 0];
  ddz_stamp = [0 -1 0; 1 0 0; 0 0 0];
  C = kron(ddx_stamp, ddx) + kron(ddy_stamp, ddy) + kron(ddz_stamp, ddz);
end
