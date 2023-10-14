function d2 = laplacian_1D(i,idx_step,idx_offset,X,N,dz)
  % e.g. find dzdz E_i with 3 poles
  % X = [E_1, dtE_1, P_11, dtP_11, P_21, dtP_21, P_31, dtP_31, E_2, dtE_2, P_12, dtP_12 ...]
  %      1    2      3     4       5     6       7     8       9    10     11    12
  % idx_step = 2*(1 + 3) = 8;
  % idx_offset = 0;
  % dzdz E_1 = (E_2 - 2*E_1)/dz^2
  %          = (X(idx_step*(i-1+1)+1) - 2*X(idx_step*(i-1)+1))/dz^2 (i = 1)
  % dzdz E_2 = (E_3 + E_1 - 2*E_2)/dz^2 
  %          = (X(idx_step*(i-1+1)+1) + X(idx_step*(i-1-1)+1) - 2*X(idx_step*(i-1)+1))/dz^2 (i = 2)
  if (i == 1)
    return (X(idx_step*(i)+1)                       - 2*X(idx_step*(i-1)+1))/dz^2;
  elseif (i == N)
    return (                    X(idx_step*(i-2)+1) - 2*X(idx_step*(i-1)+1))/dz^2;
  else
    return (X(idx_step*(i)+1) + X(idx_step*(i-2)+1) - 2*X(idx_step*(i-1)+1))/dz^2;
  end
end
