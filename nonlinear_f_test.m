clear all
close all

params = nonlinear_params();

%figure(3);
%plot_permittivity(params);
%return;

% source 
dt = params.dz/3e8; %50e-18; % s
tf = 50e-15; % s
tsteps = round(tf/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s
J = zeros(params.N,tsteps);
dtJ = zeros(params.N,tsteps);

omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
t0_J = 2*2*pi/omega_J;
dtJ(round(params.N/2),:) = ricker(5e7/params.dz*omega_J, omega_J, tsteps, dt, t0_J);
%figure(4);
%plot(tspan, dtJ(round(params.N/2),:), '-o');
%return;

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.dz,params.N); % m

E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz,1), params.N);
dtP0 = zeros(size(params.Lorentz,1), params.N);

X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

options = odeset('RelTol',1e-3,'AbsTol',1e-12);

eval_f = @(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params);
[t,X] = ode45(eval_f, tspan, X0, options);

dxdt_end = nonlinear_f(X(end,:)', 0, params);
X_end = X(end,:);
[dtE_end, dtdtE_end, dtP_end, dtdtP_end] = nonlinear_split_X(dxdt_end, params);
[E_end, ~, P_end, ~] = nonlinear_split_X(X_end, params);

gen_video = true;
plot_all = false;
if gen_video
  f = figure(1);
  set(f, 'resize', 'off', 'Position', [100 100 800 300]);
  tiledlayout(2,7);
  ax_ED = nexttile(1, [1 4]);
  ax_dtEdtD = nexttile(8, [1 4]);
  ax_J = nexttile(5, [2 3]);
  
  video = VideoWriter('nonlinear_f_test.avi'); %Create a video object
  open(video); % Open video source - restricts the use of video for your program
  
  i_list = 1:round(size(tspan,2)/200):size(tspan,2);
  E_scale = 1/1e9; % V/nm
  D_scale = 1e12/1e12; % pC/um^2
  w = round(2.5*2*pi*3e8/(omega_J*3)/params.dz);
  [E, D, dtE, dtD, P, dtP] = nonlinear_u(X(round(t0_J/dt),:)',params);
  m_E = 1.1*max(abs(E))*E_scale;
  m_D = 1.1*max(abs(D))*D_scale;
  if (m_E == 0); m_E = 1; end
  if (m_D == 0); m_D = 1; end
  for i=i_list
    [E, D, dtE, dtD, P, dtP] = nonlinear_u(X(i,:)',params);
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


    eps_J = 1e-4;
    J_f = JacobianCalculation(@(X) nonlinear_f(X,dtJ(:,i),params), X(i,:)', eps_J, size(X,2));
    set(gcf, 'CurrentAxes', ax_J);
    spy(J_f > 0, 'r');
    hold on;
    spy(J_f < 0, 'b');
    title(ax_J, sprintf("Jacobian (eps = %1.1e, x_{order} = %i)", eps_J, params.x_order));

    drawnow;
    vidFrame = getframe(gcf);
    %if i < i_list(end);
    %  clf;
    %end
    writeVideo(video,vidFrame);
  end
  
  close(video);
else
  figure(1);
  yyaxis left; plot(dtdtP_end); yyaxis right; plot(dtP_end);
end
