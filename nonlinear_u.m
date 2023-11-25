function u = nonlinear_u(t,params)
  u = zeros(params.N,1);
  if params.source_type == "ricker"
    u(round(params.N/2)) = ricker(t, params.ampl_J, params.omega_J, params.t0_J);
  elseif params.source_type == "sinusoid"
    u(round(params.N/2)) = params.ampl_J*sin(params.omega_J*t);
  end
end
