function dxdt = evalf_linear_dispersive(X,J,p)
    dxdt = zeros(4*p.N, 1);
    H = X(1:p.N);
    E = X(p.N+1:2*p.N);
    P = X(2*p.N+1:3*p.N);
    DP = X(3*p.N+1:4*p.N);
    dHdt = zeros(p.N,1);
    dEdt = zeros(p.N,1);
    dPdt = zeros(p.N,1);
    dDPdt = zeros(p.N,1);
    for i = 1:p.N
        % eqn 1: d_t H_y = 1/mu_0 d_x E_z
        if i == 1
            dHdt(i) = 1/p.mu_0*(E(i+1))/(2*p.delta_x);
        elseif i == p.N
            dHdt(i) = 1/p.mu_0*(-E(i-1))/(2*p.delta_x);
        else
            dHdt(i) = 1/p.mu_0*(E(i+1)-E(i-1))/(2*p.delta_x);
        end
        % eqn 2: d_t E_z = 1/eps_0 (d_x H_y - DP_z - J_z)
        if i == 1 
           dEdt(i) = 1/p.eps_0*((H(i+1))/(2*p.delta_x) - DP(i) - J(i));
        elseif i == p.N 
           dEdt(i) = 1/p.eps_0*(-(H(i-1))/(2*p.delta_x) - DP(i) - J(i));
        else
           dEdt(i) = 1/p.eps_0*((H(i+1)-H(i-1))/(2*p.delta_x) - DP(i) - J(i));
        end
        % eqn 3: d_t P_z = DP_z
        dPdt(i) = DP(i);
        % eqn 4: d_t DP_z = omega_chi_1^2(eps_0 chi_1 E_z - P_z - delta_chi_1/omega_chi_1^2 DP_z)
        dDPdt(i) = p.omega_chi_1^2*(p.eps_0*p.chi_1*E(i)-P(i)-(p.delta_chi_1/p.omega_chi_1^2)*DP(i));
    end
    dxdt = [dHdt; dEdt; dPdt; dDPdt];
end
