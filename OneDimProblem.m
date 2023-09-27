function [f, M] = OneDimProblem(H, E, E_dt, P, P_dt, J, N,  delta_x, eps_0, mu_0, chi_1, chi_2, omega_chi_1, delta_chi_1)
    N_equation = 5; 
    f = zeros(N_equation*N, 1);
    M = diag([ones(2*N, 1); zeros(N, 1); ones(2*N, 1)]);
    % H
    for i = 1:N
        if i == 1
            f(i) = 1/mu_0*(E(i+1))/(2*delta_x);
        elseif i == N
            f(i) = 1/mu_0*(-E(i-1))/(2*delta_x);
        else
            f(i) = 1/mu_0*(E(i+1)-E(i-1))/(2*delta_x);
        end
    end
    % E
    for i = N+1 : 2*N
        f(i) = E_dt(i); 
    end
    % conservation laws (E')
    for i = 2*N+1 : 3*N
        if i == 1 
           f(i) = -(H(i+1))/(2*delta_x) + eps_0*E_dt(i)+P_dt(i) + 2*eps_0*chi_2*E(i)*E_dt(i)+J(i);
        elseif i == N 
           f(i) = (H(i-1))/(2*delta_x) + eps_0*E_dt(i)+P_dt(i) + 2*eps_0*chi_2*E(i)*E_dt(i)+J(i);
        else
        f(i) = -(H(i+1)-H(i-1))/(2*delta_x) + eps_0*E_dt(i)+P_dt(i) + 2*eps_0*chi_2*E(i)*E_dt(i)+J(i);
        end
    end
    % P 
    for i = 3*N+1 : 4*N
        f(i) = P_dt(i);
    end 
    % P' 
    for i = 4*N+1 : 5*N
        f(i) = eps_0*omega_chi_1^2*chi_1*E(i)-omega_chi_1^2*P(i)-delta_chi_1*P_dt(i);
    end
end