function [x,converged,err_f_k,err_dx_k,err_rel_k,iterations,X] = newton(eval_f,eval_Jf,p,u,x0,newton_opts);

   % usage
   % newton(@f, @jf, p, u, x0, newton_opts);
   % eval_f: evaluates system (x, p, u)
   % eval_Jf: jacobian of system (x, p, u)
   % p: parameters
   % u: input
   % x0: initial guess
   % newton_opts: struct with options for Newton:
   %  err_f: termination condition for error on |f(x)|
   %  err_dx: termination condition on |dx|
   %  err_rel: termination condition on |dx./x|
   %  max_iter: maximum number of iterations

   k = 1;
   X(:,1) = x0;
   f = eval_f(X(:,k),p,u);
   err_f_k      = zeros(1);
   err_dx_k     = zeros(1);
   err_rel_k    = zeros(1);
   err_f_k(1)   = Inf;
   err_dx_k(1)  = Inf;
   err_rel_k(1)  = Inf;
   while k<=newton_opts.max_iter & (err_f_k(end)>newton_opts.err_f | err_dx_k(end)>newton_opts.err_dx),
      f = eval_f(X(:,k),p,u);
      Jf = eval_Jf(X(:,k),p,u);
      dx = Jf\(-f);
      X(:,k+1) = X(:,k) + dx;
      err_f_k(k+1) = norm(f,inf);
      err_dx_k(k+1) = norm(dx,inf);
      err_rel_k(k+1) = norm(dx./X(:,k),inf);
      k = k+1;
   end
   x = X(:,end);
   iterations = k-1; 
   if err_f_k(end)<=newton_opts.err_f & err_dx_k(end)<=newton_opts.err_dx & err_rel_k(end)<=newton_opts.err_rel
      converged = 1;
   else
      converged=0;
      fprintf(1, 'Newton did NOT converge! Maximum Number of Iterations reached\n');
   end
end
