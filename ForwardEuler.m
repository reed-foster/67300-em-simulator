function [X,t] = ForwardEuler(eval_f,x_start,J,p,t_start,t_stop,timestep,visualize)
% uses Forward Euler to simulate states model dx/dt=f(x,p,u)
% startin from state vector x_start at time t_start
% until time t_stop, with time intervals timestep
% eval_f is a text string defining the name of the function that evaluates f(x,p,u)
% visualize ~= 0 is an optional parameter triggering the generation of intermediate plots of the state
% 

X(:,1) = x_start;
t(1)   = t_start;

if visualize
   VisualizeState(t,X,1,'.b');
end

for n = 1 : ceil((t_stop-t_start)/timestep)
   dt       = min(timestep, (t_stop-t(n)));
   t(n+1)   = t(n) + dt;
   f        = feval(eval_f, X(:,n), J, p);
   X(:,n+1) = X(:,n) +  dt * f;
   if visualize
      VisualizeState(t,X,n+1,'.b');
   end
end
