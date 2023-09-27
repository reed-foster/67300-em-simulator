function dxdt = evalf_linear_dispersive(X,J,p)
    dxdt = zeros(4*p.N, 1);
    for i = 1:p.N
        % H <- X(1:N)
        % E <- X(p.N+1:2*p.N)
        % P <- X(2*p.N+1:3*p.N)
        % DP <- X(3*p.N+1:4*p.N)
        % eqn 1: d_t H_y = 1/mu_0 d_x E_z
        if i == 1
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N))/(2*p.delta_x);
        elseif i == p.N
            dxdt(i) = 1/p.mu_0*(-X(i-1+p.N))/(2*p.delta_x);
        else
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N)-X(i-1+p.N))/(2*p.delta_x);
        end
        % eqn 2: d_t E_z = 1/eps_0 (d_x H_y - DP_z - J_z)
        if i == 1 
           dxdt(i+p.N) = 1/p.eps_0*((X(i+1))/(2*p.delta_x) + X(i+3*p.N));
        elseif i == p.N 
           dxdt(i+p.N) = 1/p.eps_0*(-(X(i-1))/(2*p.delta_x) + X(i+3*p.N));
        else
           dxdt(i+p.N) = 1/p.eps_0*((X(i+1)-X(i-1))/(2*p.delta_x) + X(i+3*p.N));
        end
        % eqn 3: d_t P_z = DP_z
        dxdt(i+2*p.N) = X(i+3*p.N);
        % eqn 4: d_t DP_z = a E_z - b P_z - c DP_z
        dxdt(i+3*p.N) = p.omega_chi_1^2*(p.eps_0*p.chi_1*X(i+p.N)-X(i+2*p.N)-(p.delta_chi_1/p.omega_chi_1^2)*X(i+3*p.N));
    end
end
