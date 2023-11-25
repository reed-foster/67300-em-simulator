function [t,X,max_k] = trapezoid(eval_f,eval_Jf,p,u,x0,tf,dt,trap_opts,newton_opts);

   % usage
   % trapezoid(@f, @jf, p, @u, x0, tf, dt, trap_opts);
   % eval_f: evaluates system (x, p, u)
   % eval_Jf: jacobian of system (x, p, u)
   % p: parameters
   % u: input function
   % x0: initial state
   % tf: simulation stop time
   % dt: timestep
   % trap_opts: struct with options for trapezoidal integrator
   %  save_intermediate: if true, save all timesteps, otherwise only save the final state
   % newton_opts: struct with options for Newton:
   %  err_f: termination condition for error on |f(x)|
   %  err_dx: termination condition on |dx|
   %  max_iter: maximum number of iterations
   t(1) = 0;
   X(:,1) = x0;
   max_k = 1;

   for l = 1:ceil(tf/dt)
      dt_l = min(dt, tf-t(l));
      u_t = u(t(l));
      % part of f_trap that doesn't depend on x, so precompute it to save time
      gamma = X(:,l) + dt_l/2*eval_f(X(:,l),p,u_t);
      % trap function and jacobian
      f_trap = @(x,p,u) x - dt_l/2*eval_f(x,p,u) - gamma;
      Jf_trap = @(x,p,u) eye(length(x0)) - dt_l/2*eval_Jf(x,p,u);
      % call newton to solve f_trap
      [x,~,~,~,k,~] = newton(f_trap, Jf_trap, p, u_t, X(:,l), newton_opts);
      if k > max_k
         max_k = k;
      end
      % use number of iterations to dynamically adjust timestep
      X(:,l+1) = x;
      t(l+1) = t(l) + dt_l;
   end
end
