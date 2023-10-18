function d2 = laplacian(X, h);
  Xm1 = [X(2:end); 0]; 
  Xp1 = [0; X(1:end-1)];
  d2 = (Xm1 + Xp1 - 2*X)./h^2;
end
