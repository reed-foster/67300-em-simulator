function dxdt = evalf_freespace(X,J,p)
    dxdt = zeros(2*p.N, 1);
    for i = 1:p.N
        % dHdt <- dXdt(1:N)
        % dEdt = dXdt(p.N+1:2*p.N)
        % eqn 1
        if i == 1
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N))/(2*p.delta_x);
        elseif i == p.N
            dxdt(i) = 1/p.mu_0*(-X(i-1+p.N))/(2*p.delta_x);
        else
            dxdt(i) = 1/p.mu_0*(X(i+1+p.N)-X(i-1+p.N))/(2*p.delta_x);
        end
        % eqn 2
        if i == 1 
           dxdt(i+p.N) = 1/p.eps_0*((X(i+1))/(2*p.delta_x) + J(i));
        elseif i == p.N 
           dxdt(i+p.N) = 1/p.eps_0*(-(X(i-1))/(2*p.delta_x) + J(i));
        else
            dxdt(i+p.N) = 1/p.eps_0*((X(i+1)-X(i-1))/(2*p.delta_x) + J(i));
        end
    end
end