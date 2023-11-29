function [t,X,max_k] = trapezoid(eval_f,p,u,x0,tf,dt,trap_opts,newton_opts);
   % trapezoidal integrator to simulate state evolution of model dx/dt=f(x,p,u)
   %
   % usage:
   % trapezoid(@f, @jf, p, @u, x0, tf, dt, trap_opts);
   %
   % eval_f: handle to function that evaluates system f(x, p, u)
   % p: parameters
   % u: handle to input function u(t, p)
   % x0: initial state
   % tf: simulation stop time
   % dt: timestep
   % trap_opts: struct with options for trapezoidal integrator
   %  visualize: if Inf, don't visualize, otherwise, show every 'visualize'th frame
   %  save_intermediate: if true, save all results, otherwise only save final result
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
      if trap_opts.save_intermediate
         lcurr = l;
         lnext = l+1;
      else
         lcurr = 1;
         lnext = 1;
      end
      dt_l = min(dt, tf-t(lcurr));
      u_t = u(t(lcurr), p);
      % part of f_trap that doesn't depend on x, so precompute it to save time
      gamma = X(:,lcurr) + dt_l/2*eval_f(X(:,lcurr),p,u_t);
      % trap function and jacobian
      f_trap = @(x,p,u) x - dt_l/2*eval_f(x,p,u) - gamma;
      % call newton to solve f_trap
      [x,converged,err_f_k,err_dx_k,err_rel_k,k,~] = newton(f_trap, p, u_t, X(:,lcurr) + dt_l*eval_f(X(:,lcurr),p,u_t), newton_opts);
      if k > max_k
         max_k = k;
      end
      if mod(l, trap_opts.visualize) == 0
         visualize_state(X(:,lcurr),t(lcurr),p,visualize_struct);
      end
      % use number of iterations to dynamically adjust timestep
      X(:,lnext) = x;
      t(lnext) = t(lcurr) + dt_l;
   end
end
