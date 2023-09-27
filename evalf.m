function dxdt = evalf(X,J,p)
    dxdt = zeros(5*p.N, 1);
    for i = 1:p.N
        % dHdt <- dXdt(1:N)
        % dEdt = dXdt(p.N+1:2*p.N)
        % dE'dt = dXdt(2*p.N+1:3*p.N)
        % dPdt -> dXdt(3*p.N+1:4*p.N)
        % dP'dt -> dXdt(4*p.N+1:5*p.N)
        % eqn 1
        if i == 1
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N))/(2*p.delta_x);
        elseif i == p.N
            dxdt(i) = 1/p.mu_0*(-X(i-1+p.N))/(2*p.delta_x);
        else
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N)-X(i-1+p.N))/(2*p.delta_x);
        end
        % eqn 2
        dxdt(i+p.N) = X(i+2*p.N);
        % eqn 3
        if i == 1 
           dxdt(i+2*p.N) = -(X(i+1))/(2*p.delta_x) + p.eps_0*X(i+2*p.N)+X(i+4*p.N) + 2*p.eps_0*p.chi_2*X(i+p.N)*X(i+2*p.N)+J(i);
        elseif i == p.N 
           dxdt(i+2*p.N) = (X(i-1))/(2*p.delta_x) + p.eps_0*X(i+2*p.N)+X(i+4*p.N) + 2*p.eps_0*p.chi_2*X(i+p.N)*X(i+2*p.N)+J(i);
        else
            dxdt(i+2*p.N) = -(X(i+1)-X(i-1))/(2*p.delta_x) + p.eps_0*X(i+2*p.N)+X(i+4*p.N) + 2*p.eps_0*p.chi_2*X(i+p.N)*X(i+2*p.N)+J(i);
        end
        % eqn 4
        dxdt(i+3*p.N) = X(i+4*p.N);
        % eqn 5
        dxdt(i+4*p.N) = p.eps_0*p.omega_chi_1^2*p.chi_1*X(i+p.N)-p.omega_chi_1^2*X(i+3*p.N)-p.delta_chi_1*X(i+4*p.N);
    end
end
