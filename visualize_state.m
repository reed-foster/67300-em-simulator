function visualize_state(X, t, params);
  f = figure(1);
  set(f, 'resize', 'off', 'Position', [100 100 1380 820]);
  tiledlayout(2,1);
  ax_EP = nexttile(1);
  ax_dtEdtP = nexttile(2);

  xvec = linspace(0,(params.N-1)*params.dz,params.N);

  [E, dtE, P, dtP] = split_X(X./(params.X_scale), params);
  yyaxis(ax_EP, 'left');
  plot(ax_EP, xvec*1e6, E, '-o');
  ylabel(ax_EP, "field [V/m]");
  yyaxis(ax_EP, 'right');
  plot(ax_EP, xvec*1e6, sum(P,1), '-o');
  ylabel(ax_EP, "polarization density [C/m^2]");
  legend(ax_EP, "E_x(z,t)", "P_x(z,t)");
  title(ax_EP, sprintf("t = %0.3f [fs]", t*1e15));
  
  % plot time derivatives of fields
  yyaxis(ax_dtEdtP, 'left');
  plot(ax_dtEdtP, xvec*1e6, dtE, '-o');
  ylabel(ax_dtEdtP, "field [V/m/s]");
  yyaxis(ax_dtEdtP, 'right');
  plot(ax_dtEdtP, xvec*1e6, sum(dtP,1), '-o');
  ylabel(ax_dtEdtP, "polarization density [C/m^2/s]");
  legend(ax_dtEdtP, "\partial_t E_x(z,t)", "\partial_t P_x(z,t)");
  xlabel(ax_dtEdtP, "x [um]");

  drawnow;
end
