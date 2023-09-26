function C = curl_matrix_sparse(L,M,N)
  %
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
