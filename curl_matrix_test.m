close all;
L = 11;
M = 16;
N = 13;

[x,y,z] = meshgrid(-(L-1):2:(L-1),-(M-1):2:(M-1),-(N-1):2:(N-1));
Fx = 200 - (2*sin(pi*x.^2/(50)) + 2*y.^2 + 10*exp(-z.^2/100));
Fy = 200 - (10*sin(pi*x.^2/(50)) + 0.1*y.^2 + 2*exp(-z.^2/10));
Fz = 200 - (0.1*sin(pi*x.^2/(50)) + 4*y.^2 + 0.1*exp(-z.^2/800));

% matlab implementation of divergence (does a central difference in the middle and left/right-handed difference on the boundaries)
[mat_curl_x,mat_curl_y,mat_curl_z,~] = curl(x,y,z,Fx,Fy,Fz);
mat_curl_x = permute(mat_curl_x, [2 1 3]);
mat_curl_y = permute(mat_curl_y, [2 1 3]);
mat_curl_z = permute(mat_curl_z, [2 1 3]);

c = sparse(curl_matrix(L,M,N)) / 4;
raw = c*[reshape(permute(Fx,[2 1 3]),[L*M*N 1]); reshape(permute(Fy,[2 1 3]),[L*M*N 1]); reshape(permute(Fz,[2 1 3]),[L*M*N 1])];

curl_x = reshape(raw(1:L*M*N), [L M N]);
curl_y = reshape(raw(L*M*N+1:2*L*M*N), [L M N]);
curl_z = reshape(raw(2*L*M*N+1:3*L*M*N), [L M N]);
% trim boundaries, since edge calculation is slightly different from matlab's implementation
% we don't really care too much about exactly what the curl calculation will be at the boundary, since the field should be close to zero when we use a PML
% if we add periodic boundary conditions, then we can change the matrix so that the difference terms which wrap around a boundary are included.
trim = 2;

mat_curl_x_trim = mat_curl_x(trim:end-trim,trim:end-trim,trim:end-trim);
mat_curl_y_trim = mat_curl_y(trim:end-trim,trim:end-trim,trim:end-trim);
mat_curl_z_trim = mat_curl_z(trim:end-trim,trim:end-trim,trim:end-trim);
curl_x_trim = curl_x(trim:end-trim,trim:end-trim,trim:end-trim);
curl_y_trim = curl_y(trim:end-trim,trim:end-trim,trim:end-trim);
curl_z_trim = curl_z(trim:end-trim,trim:end-trim,trim:end-trim);

err_x = norm(curl_x_trim - mat_curl_x_trim, "fro");
err_y = norm(curl_y_trim - mat_curl_y_trim, "fro");
err_z = norm(curl_z_trim - mat_curl_z_trim, "fro");

disp(sprintf("Curl_x error (loop) = %d", err_x));
disp(sprintf("Curl_y error (loop) = %d", err_y));
disp(sprintf("Curl_z error (loop) = %d", err_z));

% test sparse implementation
c_sparse = curl_matrix_sparse(L,M,N) / 4;
raw_sparse = c*[reshape(permute(Fx,[2 1 3]),[L*M*N 1]); reshape(permute(Fy,[2 1 3]),[L*M*N 1]); reshape(permute(Fz,[2 1 3]),[L*M*N 1])];

curl_sparse_x = reshape(raw_sparse(1:L*M*N), [L M N]);
curl_sparse_y = reshape(raw_sparse(L*M*N+1:2*L*M*N), [L M N]);
curl_sparse_z = reshape(raw_sparse(2*L*M*N+1:3*L*M*N), [L M N]);

curl_sparse_x_trim = curl_sparse_x(trim:end-trim,trim:end-trim,trim:end-trim);
curl_sparse_y_trim = curl_sparse_y(trim:end-trim,trim:end-trim,trim:end-trim);
curl_sparse_z_trim = curl_sparse_z(trim:end-trim,trim:end-trim,trim:end-trim);

err_sparse_x = norm(curl_sparse_x_trim - mat_curl_x_trim, "fro");
err_sparse_y = norm(curl_sparse_y_trim - mat_curl_y_trim, "fro");
err_sparse_z = norm(curl_sparse_z_trim - mat_curl_z_trim, "fro");

disp(sprintf("Curl_x error (sparse) = %d", err_sparse_x));
disp(sprintf("Curl_y error (sparse) = %d", err_sparse_y));
disp(sprintf("Curl_z error (sparse) = %d", err_sparse_z));


% spy curl matrix
color_spy(curl_matrix(5,3,4));
yline(61);
yline(121);
xline(61);
xline(121);
