function u = eval_u(t,params)
  u = zeros(params.N,1);
  if params.source_type == "ricker"
    % turn on slowly
    scale = heaviside(t-params.t0_J/10) + heaviside(params.t0_J/10 - t)/(1+(1/(t/(params.t0_J/10))-1)^2);
    u(round(params.N/2)) = scale*ricker(t, params.ampl_J, params.omega_J, params.t0_J);
  elseif params.source_type == "sinusoid"
    % turn on slowly over two periods
    scale = heaviside(t-4*pi/params.omega_J) + heaviside(4*pi/params.omega_J - t)/(1+(1/(t*params.omega_J/(4*pi))-1)^2);
    u(round(params.N/2)) = scale*params.ampl_J*cos(params.omega_J*t);
  end
end
