function visualize_struct = visualize_state(X, t, params, visualize_struct);
  if visualize_struct.init
    disp("initializing visualizer");
    visualize_struct.f = figure(1);
    set(visualize_struct.f, 'resize', 'off', 'Position', [100 100 1380 820]);
    tiledlayout(2,1);
    visualize_struct.ax_E = nexttile(1);
    visualize_struct.ax_dtE = nexttile(2);
    visualize_struct.init = false;
    visualize_struct.zvec = linspace(0,(params.N-1)*params.dz,params.N);

    visualize_struct.E = plot(visualize_struct.ax_E, visualize_struct.zvec*1e6, zeros(params.N,1), '-o');
    legend(visualize_struct.ax_E, "E_x(z,t)");
    ylabel(visualize_struct.ax_E, "field [V/m]");
    xlim(visualize_struct.ax_E, [0 (params.N-1)*params.dz*1e6]);
    visualize_struct.ax_E.FontSize = 18;

    visualize_struct.dtE = plot(visualize_struct.ax_dtE, visualize_struct.zvec*1e6, zeros(params.N,1), '-o');
    legend(visualize_struct.ax_dtE, "\partial_t E_x(z,t)");
    ylabel(visualize_struct.ax_dtE, "field [V/m/s]");
    xlim(visualize_struct.ax_dtE, [0 (params.N-1)*params.dz*1e6]);
    visualize_struct.ax_dtE.FontSize = 18;

    xlabel(visualize_struct.ax_dtE, "z [um]");

  else
    figure(visualize_struct.f);
  end


  [E, dtE, ~, ~] = split_X(X./(params.X_scale), params);
  m_E = max(1.1*max(abs(E)), eps);
  m_dtE = max(1.1*max(abs(dtE)), eps);

  visualize_struct.E.YData =  E;
  ylim(visualize_struct.ax_E, [-m_E m_E]);
  
  % plot time derivatives of fields
  visualize_struct.dtE.YData = dtE;
  ylim(visualize_struct.ax_dtE, [-m_dtE m_dtE]);
  
  title(visualize_struct.ax_E, sprintf("t = %0.3f [fs]", t*1e15));

  drawnow;
end
