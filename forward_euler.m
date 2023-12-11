function [t,X] = forward_euler(eval_f,p,u,x0,tf,dt,opts)
% uses Forward Euler to simulate states model dx/dt=f(x,p,u)
%
% usage:
% forward_euler(@f, p, @u, x0, tf, dt, opts);
% 
% eval_f: a handle of the function that evaluates f(x,p,u)
% p: parameters
% u: a handle to a function u(t,p) for the time-dependent source
% x0: initial state
% tf: simulation stop time
% dt: timestep
% opts: struct with options for forward euler
%  visualize_dt: if Inf, don't visualize, otherwise show state every visualize_dt step
%  save_intermediate: if true, save all results, otherwise only save final result

X(:,1) = x0;
t(1)   = 0;

if (opts.visualize_dt < Inf)
   visualize_struct.init = true;
   visualize_struct = visualize_state(X(:,1), t(1), p, visualize_struct);
end

t_last_visualized = 0;
for n = 1 : ceil(tf/dt)
  if opts.save_intermediate
    l = n;
    lp = n+1;
  else
    l = 1;
    lp = 1;
  end
  dt_l     = min(dt, tf-t(l));
  u_t      = u(t(l), p);
  f        = eval_f(X(:,l), p, u_t);
  t(lp)   = t(l) + dt_l;
  X(:,lp) = X(:,l) +  dt_l * f;
  % visualize
  if t_last_visualized + opts.visualize_dt < t(end)
    visualize_state(X(:,end),t(end),p,visualize_struct);
    t_last_visualized = t(end);
  end
end
