function [t,X,max_k] = trapezoid(eval_f,eval_Jf,p,u,x0,tf,dt,trap_opts,newton_opts);
   % trapezoidal integrator to simulate state evolution of model dx/dt=f(x,p,u)
   %
   % usage:
   % trapezoid(@f, @jf, p, @u, x0, tf, dt, trap_opts);
   %
   % eval_f: handle to function that evaluates system f(x, p, u)
   % eval_Jf: handle to jacobian of system Jf(x, p, u)
   % p: parameters
   % u: handle to input function u(t, p)
   % x0: initial state
   % tf: simulation stop time
   % dt: timestep
   % trap_opts: struct with options for trapezoidal integrator
   %  visualize: if Inf, don't visualize, otherwise, show every 'visualize'th frame
   % newton_opts: struct with options for Newton:
   %  err_f: termination condition for error on |f(x)|
   %  err_dx: termination condition on |dx|
   %  max_iter: maximum number of iterations
   t(1) = 0;
   X(:,1) = x0;
   max_k = 1;

   if (trap_opts.visualize < Inf)
      visualize_struct.init = true;
      visualize_struct = visualize_state(X(:,1), t(1), p, visualize_struct);
   end

   for l = 1:ceil(tf/dt)
      dt_l = min(dt, tf-t(l));
      u_t = u(t(l), p);
      % part of f_trap that doesn't depend on x, so precompute it to save time
      gamma = X(:,l) + dt_l/2*eval_f(X(:,l),p,u_t);
      % trap function and jacobian
      f_trap = @(x,p,u) x - dt_l/2*eval_f(x,p,u) - gamma;
      Jf_trap = @(x,p,u) eye(length(x0)) - dt_l/2*eval_Jf(x,p,u);
      % call newton to solve f_trap
      [x,converged,err_f_k,err_dx_k,err_rel_k,k,~] = newton(f_trap, Jf_trap, p, u_t, X(:,l) + dt_l*eval_f(X(:,l),p,u_t), newton_opts);
      if k > max_k
         max_k = k;
      end
      if mod(l, trap_opts.visualize) == 0
         visualize_state(X(:,l),t(l),p,visualize_struct);
      end
      % use number of iterations to dynamically adjust timestep
      X(:,l+1) = x;
      t(l+1) = t(l) + dt_l;
   end
end
