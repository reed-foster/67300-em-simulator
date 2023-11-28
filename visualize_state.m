function visualize_struct = visualize_state(X, t, params, visualize_struct);
  if visualize_struct.init
    disp("initializing visualizer");
    visualize_struct.f = figure(1);
    set(visualize_struct.f, 'resize', 'off', 'Position', [100 100 1380 820]);
    tiledlayout(2,1);
    visualize_struct.ax_EP = nexttile(1);
    visualize_struct.ax_dtEdtP = nexttile(2);
    visualize_struct.init = false;
    visualize_struct.zvec = linspace(0,(params.N-1)*params.dz,params.N);
  else
    figure(visualize_struct.f);
  end


  [E, dtE, P, dtP] = split_X(X./(params.X_scale), params);
  P_tot = sum(P,1);
  dtP_tot = sum(dtP,1);
  m_E = max(1.1*max(abs(E)), eps);
  m_dtE = max(1.1*max(abs(dtE)), eps);
  m_P = max(1.1*max(abs(P_tot)), eps);
  m_dtP = max(1.1*max(abs(dtP_tot)), eps);

  yyaxis(visualize_struct.ax_EP, 'left');
  plot(visualize_struct.ax_EP, visualize_struct.zvec*1e6, E, '-o');
  ylabel(visualize_struct.ax_EP, "field [V/m]");
  ylim(visualize_struct.ax_EP, [-m_E m_E]);

  yyaxis(visualize_struct.ax_EP, 'right');
  plot(visualize_struct.ax_EP, visualize_struct.zvec*1e6, P_tot, '-o');
  ylabel(visualize_struct.ax_EP, "polarization density [C/m^2]");
  ylim(visualize_struct.ax_EP, [-m_P m_P]);

  legend(visualize_struct.ax_EP, "E_x(z,t)", "P_x(z,t)");
  title(visualize_struct.ax_EP, sprintf("t = %0.3f [fs]", t*1e15));
  
  % plot time derivatives of fields
  yyaxis(visualize_struct.ax_dtEdtP, 'left');
  plot(visualize_struct.ax_dtEdtP, visualize_struct.zvec*1e6, dtE, '-o');
  ylabel(visualize_struct.ax_dtEdtP, "field [V/m/s]");
  ylim(visualize_struct.ax_dtEdtP, [-m_dtE m_dtE]);

  yyaxis(visualize_struct.ax_dtEdtP, 'right');
  plot(visualize_struct.ax_dtEdtP, visualize_struct.zvec*1e6, dtP_tot, '-o');
  ylabel(visualize_struct.ax_dtEdtP, "polarization density [C/m^2/s]");
  ylim(visualize_struct.ax_dtEdtP, [-m_dtP m_dtP]);

  legend(visualize_struct.ax_dtEdtP, "\partial_t E_x(z,t)", "\partial_t P_x(z,t)");
  xlabel(visualize_struct.ax_dtEdtP, "z [um]");

  drawnow;
end
