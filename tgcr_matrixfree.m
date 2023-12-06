function [x, r_norms, k] = tgcr_matrixfree(eval_f,xf,b,opts);
  % generalized conjugate residual method for solving [df/dx] x = b
  % using matrix-free technique
  % usage
  % tgcr_matrixfree(@f, xf, b, opts)
  % eval_f: evaluates system at x
  % xf: state to evaluate Jacobian [df/dx]
  % b: right hand side of linear system to be solved
  % opts: struct with options for GCR
  %  err_b: termination condition on error norm(b - Ax)/norm(b)
  %  max_iter: maximum number of iterations
  %  eps_x: finite difference perturbation for matrix-free directional derivative

  % initial guess for x is zero
  x = zeros(size(b));

  % set initial residual to b - Ax^0 = b (since x^0 = 0)
  r = b;
  r_norms(1) = norm(r,2);
  k = 1;
  
  % step size for finite difference approximation of Ap(:,k)
  % gets normalized by L2 norm: ||p(:,k)||_2
  epsilon_xf = opts.eps_x*norm(xf,inf);
  if epsilon_xf == 0
    epsilon_xf = eps;
  end

  while (r_norms(k)/r_norms(1) > opts.err_b) & (k <= opts.max_iter)
    % use the residual as the first guess for the next search direction
    % and compute its image in b-space
    p(:,k) = r;
    
    % matrix-free approximation for Ap(:,k) = A*p(:,k) (since A is a Jacobian)
    epsilon_norm = epsilon_xf / r_norms(k);
    Ap(:,k) = (eval_f(xf + p(:,k)*epsilon_norm)- eval_f(xf))/epsilon_norm;

    % orthogonalize Ap to previous Ap vectors
    % orthogonalize p vectors A^TA-orthogonal to previous p vectors
    % (since p is in x-space, we just orthogonalize in b-space)
    if k > 1
      for j = 1:k-1
        % Jacobian of our system is not symmetric, so we need to
        % orthogonalize against all Ap(:,j < k)
        beta = Ap(:,k)' * Ap(:,j); % projection of new Ap onto previous Ap
        p(:,k) = p(:,k) - beta * p(:,j); % orthogonalized in b-space
        Ap(:,k) = Ap(:,k) - beta * Ap(:,j);
      end
    end

    % make Ap unit length, and scale p so that it's unit length in b-space
    norm_Ap = norm(Ap(:,k),2);
    Ap(:,k) = Ap(:,k)/norm_Ap;
    p(:,k) = p(:,k)/norm_Ap;

    % determine optimal change in x along p by projecting r onto Ap
    alpha = r' * Ap(:,k);

    % update x and r
    x = x + alpha * p(:,k);
    r = r - alpha * Ap(:,k);

    % save the norm of r
    r_norms(k+1) = norm(r,2);

    k = k + 1;
  end
  if r_norms(k) > (opts.err_b * r_norms(1))
    fprintf(1, 'GCR did not converge! Maximum number of iterations reached\n');
    x = [];
  end
  % normalize r_norms with respect to initial residual norm
  r_norms = r_norms / r_norms(1);
end
