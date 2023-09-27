function dxdt = OneDimProblem(H, E, E_dt, P, P_dt, J, N,  delta_x, eps_0, mu_0, chi_1, chi_2, omega_chi_1, delta_chi_1)
    N_equation = 5; 
    dxdt = zeros(N_equation*N, 1);
    for i = 1:N
        % dHdt <- dXdt(1:N)
        if i == 1
            dxdt(i) = 1/mu_0*(E(i+1))/(2*delta_x);
        elseif i == N
            dxdt(i) = 1/mu_0*(-E(i-1))/(2*delta_x);
        else
            dxdt(i) = 1/mu_0*(E(i+1)-E(i-1))/(2*delta_x);
        end
        % dEdt = dXdt(N+1:2*N)
        dxdt(i+N) = E_dt(i); 
        % dE'dt = dXdt(2*N+1:3*N)
        if i == 1 
           dxdt(i+2*N) = -(H(i+1))/(2*delta_x) + eps_0*E_dt(i)+P_dt(i) + 2*eps_0*chi_2*E(i)*E_dt(i)+J(i);
        elseif i == N 
           dxdt(i+2*N) = (H(i-1))/(2*delta_x) + eps_0*E_dt(i)+P_dt(i) + 2*eps_0*chi_2*E(i)*E_dt(i)+J(i);
        else
        dxdt(i+2*N) = -(H(i+1)-H(i-1))/(2*delta_x) + eps_0*E_dt(i)+P_dt(i) + 2*eps_0*chi_2*E(i)*E_dt(i)+J(i);
        end
        % dPdt -> dXdt(3*N+1:4*N)
        dxdt(i+3*N) = P_dt(i);
        % dP'dt -> dXdt(4*N+1:5*N)
        dxdt(i+4*N) = eps_0*omega_chi_1^2*chi_1*E(i)-omega_chi_1^2*P(i)-delta_chi_1*P_dt(i);
    end
end
