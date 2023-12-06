function [x,converged,err_f_k,err_dx_k,err_rel_k,iterations,max_gcr_iterations,X] = newton(eval_f,p,u,x0,opts);

   % usage
   % newton(@f, p, u, x0, opts);
   % eval_f: evaluates system (x, p, u)
   % p: parameters
   % u: input
   % x0: initial guess
   % opts: struct with options for Newton:
   %  err_f: termination condition for error on |f(x)|
   %  err_dx: termination condition on |dx|
   %  err_rel: termination condition on |dx./x|
   %  max_iter: maximum number of iterations
   %  matrix_free: if true, use matrix-free GCR method, otherwise use finite-difference Jacobian
   %  err_gcr: termination condition for matrix-free GCR method on error norm(b - Ax)/norm(b)
   %  eps_fd: finite difference perturbation for matrix-free directional derivative and finite difference Jacobian

   gcr_opts.err_b = opts.err_gcr;
   gcr_opts.max_iter = size(x0,1) * 1.1;
   gcr_opts.eps_x = opts.eps_fd;

   k = 1;
   X(:,1) = x0;
   f = eval_f(X(:,k),p,u);
   err_f_k      = zeros(1);
   err_dx_k     = zeros(1);
   err_rel_k    = zeros(1);
   err_f_k(1)   = Inf;
   err_dx_k(1)  = Inf;
   err_rel_k(1)  = Inf;
   max_gcr_iterations = 1;
   while k<=opts.max_iter & (err_f_k(end)>opts.err_f | err_dx_k(end)>opts.err_dx | err_rel_k(end)>opts.err_rel),
      f = eval_f(X(:,k),p,u);
      if opts.matrix_free
         [dx, r_gcr, k_gcr] = tgcr_matrixfree(@(x) eval_f(x,p,u), X(:,k), -f, gcr_opts);
         if (k_gcr > max_gcr_iterations)
            max_gcr_iterations = k_gcr;
         end
      else
         eps_J = max(opts.eps_fd*norm(X(:,k),inf), eps);
         Jf = JacobianCalculation(@(x) eval_f(x,p,u), X(:,k), eps_J, size(x0,1));
         dx = Jf\(-f);
      end
      X(:,k+1) = X(:,k) + dx;
      err_f_k(k+1) = norm(f,inf);
      err_dx_k(k+1) = norm(dx,inf);
      err_rel_k(k+1) = norm(dx,inf)./(norm(X(:,k),inf) + eps);
      k = k+1;
   end
   x = X(:,end);
   iterations = k-1; 
   if err_f_k(end)<=opts.err_f & err_dx_k(end)<=opts.err_dx & err_rel_k(end)<=opts.err_rel
      converged = 1;
   else
      converged=0;
      fprintf(1, 'Newton did NOT converge! Maximum Number of Iterations reached\n');
   end
end
