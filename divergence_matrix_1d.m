function div = divergence_matrix_1d(N)
  % returns the 1D central finite-difference approximation of the divergence of a vector field
  div = diag(ones(1,N-1),1) + diag(-ones(1,N-1),-1);
end
