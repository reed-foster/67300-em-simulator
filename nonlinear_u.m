function u = nonlinear_u(t,params)
  u = zeros(params.N,1);
  u(round(params.N/2)) = ricker(t, params.ampl_J, params.omega_J, params.t0_J);
end
