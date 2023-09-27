function dxdt = evalf(X,J,p)
    dxdt = zeros(5*p.N, 1);
    for i = 1:p.N
        % H <- X(1:N)
        % E <- X(p.N+1:2*p.N)
        % DE <- X(2*p.N+1:3*p.N)
        % P <- X(3*p.N+1:4*p.N)
        % DP <- X(4*p.N+1:5*p.N)
        % eqn 1: d_t H_y = 1/mu_0 d_x E_z
        if i == 1
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N))/(2*p.delta_x);
        elseif i == p.N
            dxdt(i) = 1/p.mu_0*(-X(i-1+p.N))/(2*p.delta_x);
        else
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N)-X(i-1+p.N))/(2*p.delta_x);
        end
        % eqn 2: d_t E_z = DE_z
        dxdt(i+p.N) = X(i+2*p.N);
        % eqn 3: 0 = - d_x H_y + d_t (eps_0 E_z + P_z + eps_0 chi_2 E_z E_z) + J
        % eqn 3: 0 = - d_x H_y + d_t (eps_0 E_z + P_z)
        if i == 1 
           dxdt(i+2*p.N) = -(X(i+1))/(2*p.delta_x) + p.eps_0*X(i+2*p.N)+X(i+4*p.N);% + 2*p.eps_0*p.chi_2*X(i+p.N)*X(i+2*p.N)+J(i);
        elseif i == p.N 
           dxdt(i+2*p.N) = (X(i-1))/(2*p.delta_x) + p.eps_0*X(i+2*p.N)+X(i+4*p.N);% + 2*p.eps_0*p.chi_2*X(i+p.N)*X(i+2*p.N)+J(i);
        else
            dxdt(i+2*p.N) = -(X(i+1)-X(i-1))/(2*p.delta_x) + p.eps_0*X(i+2*p.N)+X(i+4*p.N);% + 2*p.eps_0*p.chi_2*X(i+p.N)*X(i+2*p.N)+J(i);
        end
        % eqn 4: d_t P_z = DP_z
        dxdt(i+3*p.N) = X(i+4*p.N);
        % eqn 5: d_t DP_z = a E_z - b P_z - c DP_z
        dxdt(i+4*p.N) = p.eps_0*p.omega_chi_1^2*p.chi_1*X(i+p.N)-p.omega_chi_1^2*X(i+3*p.N)-p.delta_chi_1*X(i+4*p.N);
    end
end
