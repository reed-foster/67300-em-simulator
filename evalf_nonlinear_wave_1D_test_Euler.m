close all

params = nonlinear_params();

%parameters 
params.N = 2000; % discretization
params.delta_x = 100e-9; % m

% nodal variables
E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz, 1), params.N);
dtP0 = zeros(size(params.Lorentz, 1), params.N);
X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.delta_x,params.N); % m

t_start = 0;
t_stop = 80e-15;
timestep = 10e-18;
dt = timestep;
tsteps = round(t_stop/timestep);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s

% source 
dtJ = zeros(params.N,round((t_stop-t_start)/timestep));
omega_J = 2*pi*3e8/(1.55e-6*3); % angular frequency for 1.55um
t0_J = 1.5*2*pi/omega_J;
dtJ(round(params.N/2),:) = ricker(2e8/params.delta_x*omega_J, omega_J, tsteps, timestep, t0_J);

[X,t] = ForwardEuler(@nonlinear_f, X0, dtJ, params, t_start,t_stop,timestep,0);

gen_video = true;
plot_all = false;
if gen_video
  figure(1);
  
  video = VideoWriter('evalf_nonliear_wave_1D_test_Euler.avi'); %Create a video object
  open(video); % Open video source - restricts the use of video for your program
  
  i_list = 1:round(size(tspan,2)/200):size(tspan,2);
  E_scale = 1/1e9; % V/nm
  D_scale = 1e12/1e12; % pC/um^2
  w = round(2.5*2*pi*3e8/(omega_J*3)/params.dz);
  [E, D, dtE, dtD, P, dtP] = nonlinear_u(X(:,round(t0_J/dt)),params);
  m_E = 1.1*max(abs(E))*E_scale;
  m_D = 1.1*max(abs(D))*D_scale;
  if (m_E == 0); m_E = 1; end
  if (m_D == 0); m_D = 1; end
  for i=i_list
    [E, D, dtE, dtD, P, dtP] = nonlinear_u(X(:,i),params);
    subplot(2,1,1);
    yyaxis left;
    plot(x*1e6, E*E_scale, '-o');
    ylabel("field [V/nm]");
    ylim([-m_E m_E]);
    yyaxis right;
    plot(x*1e6, D*D_scale, '-o');
    ylabel("displacement [pC/um^2]");
    ylim([-m_D m_D]);
    legend("E_x(z,t)", "D_x(z,t)");
    title(sprintf("t = %0.3f [fs]", tspan(i)*1e15));

    subplot(2,1,2);
    % zoom in on the forward-propagating pulse
    yyaxis left;
    [~,n0] = max(E(round(params.N/2):end));
    n0 = n0 + round(params.N/2) + round(w/2);
    xzoom = n0*params.dz+x(n0-w:n0+w); % m
    Ezoom = E(n0-w:n0+w);
    Dzoom = D(n0-w:n0+w);
    m_Ezoom = 1.5*max(abs(Ezoom))*E_scale;
    m_Dzoom = 1.5*max(abs(Dzoom))*D_scale;
    if (m_Ezoom == 0); m_Ezoom = 1; end
    if (m_Dzoom == 0); m_Dzoom = 1; end
    plot(xzoom*1e6, Ezoom*E_scale, '-o');
    ylabel("field [V/nm]");
    ylim([-m_Ezoom m_Ezoom]);
    xlim([xzoom(1)*1e6, xzoom(end)*1e6]);
    yyaxis right;
    plot(xzoom*1e6, Dzoom*D_scale, '-o');
    ylabel("displacement [pC/um^2]");
    ylim([-m_Dzoom m_Dzoom]);
    xlim([xzoom(1)*1e6, xzoom(end)*1e6]);
    legend("E_x(z,t)", "D_x(z,t)");
    
    xlabel("x [um]");
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
