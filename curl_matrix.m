function C = curl_matrix(L,M,N)
  % returns the central finite difference approximation of the curl of a vector field in R3
  % the matrix transform
  % the matrix transforms a 1D state vector ordered first by increasing spatial coordinate along the x,y,z axes (respectively),
  % then by vector direction e.g. V_x(0,0,0) V_x(1,0,0) ... V_x(0,1,0) ... V_y(0,0,0) ...
  % L is the number of points in space in the x direction, M the number of points in the y direction, and N in the z direction
  for i=1:L
    for j=1:M
      for k=1:N
        ijk_flat = (i-1) + L*(j-1) + L*M*(k-1) + 1;
        % curl_x
        if k < N;
          % delta_k = 1
          C(ijk_flat, L*M*N + ijk_flat + L*M) = -1;
        end
        if k > 1
          % delta_k = -1
          C(ijk_flat, L*M*N + ijk_flat - L*M) = 1;
        end
        if j < M
          % delta_j = 1
          C(ijk_flat, 2*L*M*N + ijk_flat + L) = 1;
        end
        if j > 1
          % delta_j = -1
          C(ijk_flat, 2*L*M*N + ijk_flat - L) = -1;
        end
        % curl_y
        if k < N
          % delta_k = 1
          C(ijk_flat + L*M*N, ijk_flat + L*M) = 1;
        end
        if k > 1
          % delta_k = -1
          C(ijk_flat + L*M*N, ijk_flat - L*M) = -1;
        end
        if i < L
          % delta_i = 1
          C(ijk_flat + L*M*N, 2*L*M*N + ijk_flat + 1) = -1;
        end
        if i > 1
          % delta_i = -1
          C(ijk_flat + L*M*N, 2*L*M*N + ijk_flat - 1) = 1;
        end
        % curl_z
        if j < M
          % delta_j = 1
          C(ijk_flat + 2*L*M*N, ijk_flat + L) = -1;
        end
        if j > 1
          % delta_j = -1
          C(ijk_flat + 2*L*M*N, ijk_flat - L) = 1;
        end
        if i < L
          % delta_i = 1
          C(ijk_flat + 2*L*M*N, L*M*N + ijk_flat + 1) = 1;
        end
        if i > 1
          % delta_i = -1
          C(ijk_flat + 2*L*M*N, L*M*N + ijk_flat - 1) = -1;
        end
      end
    end
  end
end
