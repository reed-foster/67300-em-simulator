clear all
close all

params = nonlinear_params();

% source 
dt = 2e-18; % s
tf = 100e-15; % s
tsteps = round(tf/dt);
tspan = linspace(0,(tsteps-1)*dt,tsteps); % s
J = zeros(params.N,tsteps);
dtJ = zeros(params.N,tsteps);

omega_J = 2*pi*6.5e13; % Hz
%dtJ_exp = 1e12*(exp(-(tspan-3/omega_J).^2/(2*(1/omega_J)^2)) - exp(-(tspan-7/omega_J).^2/(2*(1/omega_J)^2)));
dtJ_dipole = 1e3*omega_J*cos(omega_J*tspan).*exp(-(tspan-4*pi/omega_J).^2./(2*(pi/omega_J)^2));
dtJ(round(params.N/2),:) = dtJ_dipole;
%plot(tspan, dtJ(round(params.N/2),:), '-o');
%return;

% simulation setup and initial conditions
x = linspace(0,(params.N-1)*params.dz,params.N); % m

%E0 = gaussian_start(2e8, 1e-6, 5e-6, 1.55e-6, params.N, params.dz, 0);
%dtE0 = gaussian_start(2e8*1.5e17, 1e-6, 5e-6, 1.55e-6, params.N, params.dz, pi/2);
%P0 = params.eps_0*params.Lorentz(:,1)*E0;
%dtP0 = zeros(size(params.Lorentz,1), params.N);

E0 = zeros(params.N,1);
dtE0 = zeros(params.N, 1);
P0 = zeros(size(params.Lorentz,1), params.N);
dtP0 = zeros(size(params.Lorentz,1), params.N);

X0 = nonlinear_generate_X(E0, dtE0, P0, dtP0, params);


%disp_progress = @(t,y,flag)fprintf('t= %s \n',mat2str(t))*0;
%options = odeset('RelTol',1e-3,'AbsTol',1e-6,'OutputFcn',disp_progress);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);

[t,X] = ode45(@(t,X) nonlinear_f(X,dtJ(:,round(t/dt+0.5)),params), tspan, X0, options);

figure(1);

video = VideoWriter('nonlinear_f_test.avi'); %Create a video object
open(video); % Open video source - restricts the use of video for your program

i_list = 1:round(size(tspan,2)/100):size(tspan,2);
for i=i_list
  [E, D, dtE, dtD] = nonlinear_u(X(i,:)',params);
  subplot(2,1,1);
  yyaxis left;
  %plot(x*1e6, E/1e9, '-o');
  plot(x*1e6, E, '-o');
  ylabel("field [V/m]");
  %ylim([-0.2 0.2]);
  yyaxis right;
  %plot(x*1e6, D*1e3, '-o');
  %ylabel("displacement [nC/um^2]");
  plot(x*1e6, D, '-o');
  ylabel("displacement [C/m^2]");
  %ylim([-5 5]);

  subplot(2,1,2);
  yyaxis left;
  %plot(x*1e6, dtE/1e9/1e15, '-o');
  %ylabel("d/dt field [V/nm/fs]");
  plot(x*1e6, dtE, '-o');
  ylabel("d/dt field [V/m/s]");
  %ylim([-50 50]);
  yyaxis right;
  %plot(x*1e6, dtD*1e3/1e15, '-o');
  %ylabel("d/dt displacement [nC/um^2/fs]");
  plot(x*1e6, dtD, '-o');
  ylabel("d/dt displacement [C/m^2/s]");
  %ylim([-0.5 0.5]);
  
  legend("E_x(z,t)", "D_x(z,t)");
  xlabel("x [um]");
  title(sprintf("t = %0.3f [fs]", tspan(i)*1e15));
  drawnow;
  vidFrame = getframe(gcf);
  if i < i_list(end);
    clf;
  end
  writeVideo(video,vidFrame);
end

close(video);
