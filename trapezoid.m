function [t,X] = trapezoid(eval_f,p,u,x0,tf,dt_max,Jf0,trap_opts,newton_opts);
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
  % dt_max: maximum timestep. If adaptive timestepping is disabled,
  %   this is used as a fixed timestep
  % Jf0: Jacobian with zero field
  % trap_opts: struct with options for trapezoidal integrator
  %  visualize_dt: if Inf, don't visualize, otherwise, show every data every visualize_dt step
  %  save_intermediate: if true, save all results, otherwise only save final result
  %  adaptive_timestep: if true, increase the timestep adaptively without exceeding dt
  %  linear_only: if true, don't treat any nonlinearities and solve directly with a linear system solve
  % newton_opts: struct with options for Newton:
  %  err_f: termination condition for error on |f(x)|
  %  err_dx: termination condition on |dx|
  %  max_iter: maximum number of iterations
  %  preconditioner: if true, use preconditioning matrix calculated from Jf0 for GCR
  t(1) = 0;
  X(:,1) = x0;
  max_k = 1;

  if (trap_opts.visualize_dt < Inf)
    visualize_struct.init = true;
    if sum(contains(fieldnames(trap_opts), 'movie_name'))
      visualize_struct.save_movie = true;
    else
      visualize_struct.save_movie = false;
    end
    visualize_struct = visualize_state(X(:,1), t(1), p, visualize_struct);
  end

  % preconditioning matrix for GCR
  T_inv = trap_newton_gcr_preconditioner(Jf0, p);
  if trap_opts.linear_only == true
      [L,U,P] = lu(Jf0_trap);
  end

  if (trap_opts.adaptive_timestep)
    dt = dt_max/10;
  else
    dt = dt_max;
  end
  t_last_visualized = 0;
  while t(end) < tf
    dt_l = min(dt, tf-t(end));
    u_t = u(t(end), p);
    % part of f_trap that doesn't depend on x, so precompute it to save time
    gamma = X(:,end) + dt_l/2*eval_f(X(:,end),p,u_t);
    % trap function and jacobian
    f_trap = @(x,p,u) x - dt_l/2*eval_f(x,p,u) - gamma;
    if trap_opts.linear_only == true
      % Solve the system directly
      b = -f_trap(X(:,end),p,u_t);
      y = L\(P*b);
      dx = U\y;
      x = X(:,end) + dx;
    else
      % call newton to solve f_trap
      [x,converged,err_f_k,err_dx_k,err_rel_k,k,k_gcr,~] = newton(f_trap, p, u_t, X(:,end) + dt_l*eval_f(X(:,end),p,u_t), T_inv, newton_opts);
      % use number of iterations to dynamically adjust timestep
      if (trap_opts.adaptive_timestep)
        if (k_gcr > 40)
          dt = dt/5;
        elseif (k_gcr > 25)
          dt = dt/2;
        elseif (k_gcr < 20)
          dt = min(dt*1.5,dt_max);
        end
      end
      if trap_opts.print_debug
        disp(num2str([k dt], "Newton took %d iterations, setting timestep to %d"));
        disp(num2str([k_gcr dt], "GCR took max %d iterations, setting timestep to %d"));
      end
    end
    % update X
    if trap_opts.save_intermediate
      X(:,end+1) = x;
      t(end+1) = t(end) + dt_l;
    else
      X(:,end) = x;
      t(end) = t(end) + dt_l;
    end
    % visualize
    if t_last_visualized + trap_opts.visualize_dt <= t(end)
      visualize_struct = visualize_state(X(:,end),t(end),p,visualize_struct);
      t_last_visualized = t(end) - mod(t(end), trap_opts.visualize_dt);
    end
  end

  if (visualize_struct.save_movie)
    disp(num2str(length(visualize_struct.frames), "saving movie with %0d frames"));
    writerObj = VideoWriter(trap_opts.movie_name);
    open(writerObj);
    for i = 1:length(visualize_struct.frames)
      writeVideo(writerObj, visualize_struct.frames(i).cdata);
    end
    close(writerObj);
  end
end
