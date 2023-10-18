function nonlinear_f_plot(x, tspan, X, params, zoom_width, gen_video, plot_jacobian, frame_dec, fname);
  f = figure(1);
  if plot_jacobian
    set(f, 'resize', 'off', 'Position', [100 100 800 300]);
    tiledlayout(2,7);
    ax_ED = nexttile(1, [1 4]);
    ax_dtEdtD = nexttile(8, [1 4]);
    ax_J = nexttile(5, [2 3]);
  else
    set(f, 'resize', 'off', 'Position', [100 100 1380 820]);
    tiledlayout(2,1);
    ax_ED = nexttile(1);
    ax_dtEdtD = nexttile(2);
  end
  
  video = VideoWriter(fname); %Create a video object
  open(video); % Open video source - restricts the use of video for your program
  
  i_list = 1:frame_dec:size(tspan,2);
  E_scale = 1/1e9; % V/nm
  D_scale = 1e12/1e12; % pC/um^2
  w = round(zoom_width/params.dz);
  m_E = 0;
  m_D = 0;
  for i=1:size(X,2)
    [E, D, dtE, dtD, P, dtP] = nonlinear_u(X(:,i),params);
    m_E = max(1.1*max(abs(E))*E_scale, m_E);
    m_D = max(1.1*max(abs(D))*D_scale, m_D);
  end
  if (m_E == 0); m_E = 1; end
  if (m_D == 0); m_D = 1; end
  for i=i_list
    [E, D, dtE, dtD, P, dtP] = nonlinear_u(X(:,i),params);
    yyaxis(ax_ED, 'left');
    plot(ax_ED, x*1e6, E*E_scale, '-o');
    ylabel(ax_ED, "field [V/nm]");
    ylim(ax_ED, [-m_E m_E]);
    yyaxis(ax_ED, 'right');
    plot(ax_ED, x*1e6, D*D_scale, '-o');
    ylabel(ax_ED, "displacement [pC/um^2]");
    ylim(ax_ED, [-m_D m_D]);
    legend(ax_ED, "E_x(z,t)", "D_x(z,t)");
    title(ax_ED, sprintf("t = %0.3f [fs]", tspan(i)*1e15));
  
    % zoom in on the forward-propagating pulse
    yyaxis(ax_dtEdtD, 'left');
    [m,n0] = max(E(round(params.N/2):end));
    if (m < 0.001*m_E)
      n0 = round(params.N/2);
    else
      n0 = n0 + round(params.N/2) + round(w/2);
    end
    xzoom = x(max(n0-w,1):min(n0+w,params.N)); % m
    Ezoom = E(max(n0-w,1):min(n0+w,params.N));
    Dzoom = D(max(n0-w,1):min(n0+w,params.N));
    m_Ezoom = max(1.5*max(abs(Ezoom))*E_scale, 0.1*m_E);
    m_Dzoom = max(1.5*max(abs(Dzoom))*D_scale, 0.1*m_D);
    if (m_Ezoom == 0); m_Ezoom = 1; end
    if (m_Dzoom == 0); m_Dzoom = 1; end
    plot(ax_dtEdtD, xzoom*1e6, Ezoom*E_scale, '-o');
    ylabel(ax_dtEdtD, "field [V/nm]");
    ylim(ax_dtEdtD, [-m_Ezoom m_Ezoom]);
    xlim(ax_dtEdtD, [xzoom(1)*1e6, xzoom(end)*1e6]);
    yyaxis(ax_dtEdtD, 'right');
    plot(ax_dtEdtD, xzoom*1e6, Dzoom*D_scale, '-o');
    ylabel(ax_dtEdtD, "displacement [pC/um^2]");
    ylim(ax_dtEdtD, [-m_Dzoom m_Dzoom]);
    xlim(ax_dtEdtD, [xzoom(1)*1e6, xzoom(end)*1e6]);
    legend(ax_dtEdtD, "E_x(z,t)", "D_x(z,t)");
    
    xlabel(ax_dtEdtD, "x [um]");
  
    if (plot_jacobian)
      eps_J = 1e-4;
      J_f = JacobianCalculation(@(X) nonlinear_f(X,dtJ(:,i),params), X(:,i), eps_J, size(X,1));
      set(gcf, 'CurrentAxes', ax_J);
      cla(ax_J);
      spy(J_f > 0, 'r');
      hold on;
      spy(J_f < 0, 'b');
      title(ax_J, sprintf("Jacobian (eps = %1.1e, x_{order} = %i)", eps_J, params.x_order));
    end
  
    drawnow;
    vidFrame = getframe(gcf);
    writeVideo(video,vidFrame);
  end
  
  close(video);
end
