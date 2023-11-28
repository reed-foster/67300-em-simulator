clear all
close all

p = simulation_params();

% override timestep
p.dt = 1e-20;
% with Lorentizan Q = 1000:
%   p.Lorentz = [4 2*pi*1e8 2*pi*1e11; ...
%                2 2*pi*5e11 2*pi*5e14];
% excitation frequency is 193 THz, very close to resonance at 500 THz
% t_unstable = 1e-19 for forward_euler and low field amplitude (e.g. 5e-3)
% t_unstable = <1e-20 for forward_euler and high field amplitude (e.g. 5e7)
% with Lorentzian Q = 10:
%   p.Lorentz = [4 2*pi*1e10 2*pi*1e11; ...
%                2 2*pi*3e14 2*pi*3e15];
% excitation frequency is 193 THz, much further from resonance at 3 PHz
% t_unstable = 1e-19 for forward_euler at all field strengths

% initial state
E0 = zeros(p.N, 1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);

if ~isfile("reference_accuracy.mat");
  opts.visualize = Inf;
  opts.save_intermediate = false;
  dt = logspace(-15,-21,7);
  ampl = logspace(-4,8,4);
  X_ref = zeros(size(X0,1),length(ampl),length(dt));
  ref_confidence = zeros(length(ampl),length(dt));
  disp("computing reference solutions X_ref");
  for i = 1:length(ampl)
    p.ampl_J = ampl(i)/p.dz*p.omega_J;
    disp(strcat("setting dtJ_x amplitude to ", num2str(p.ampl_J), " [A/m/s]"));
    for j = 1:length(dt)
      disp(strcat("computing reference for timestep dt = ", num2str(dt(j)*1e18), " [as]"));
      tic;
      [t,X] = forward_euler(@eval_f, p, @eval_u, X0, p.tf, dt(j), opts);
      toc;
      X_ref(:,i,j) = X;
      if j == 1 % don't have a reference for first timestep
        ref_confidence(i,j) = Inf;
      else
        % calculate confidence based on relative error of E and P
        [E,~,P,~] = split_X(X, p);
        [E_ref,~,P_ref,~] = split_X(X_ref(:,i,j-1), p);
        E_err = max(abs(E - E_ref)./max(abs(E_ref)));
        P_err = max(abs(sum(P,1) - sum(P_ref,1))./max(abs(sum(P_ref,1))));
        ref_confidence(i,j) = max(E_err, P_err);
      end
      save("reference_accuracy.mat","X_ref","ref_confidence","dt","ampl","p");
    end
  end
else
  load("reference_accuracy.mat");
  % X_ref, ref_confidence, dt, ampl, p
end

for i=1:length(ampl)
  figure;
  zvec = linspace(0, (p.N-1)*p.dz, p.N);
  [E_ref, ~, P_ref, ~] = split_X(X_ref(:,i,length(dt))./p.X_scale, p);
  m_E = 1.1*max(abs(E_ref)) + eps;
  m_P = 1.1*max(abs(sum(P_ref,1))) + eps;
  for j = 1:length(dt)
    [E, ~, P, ~] = split_X(X_ref(:,i,j)./p.X_scale, p);
    subplot(2,1,1);
    plot(zvec*1e6, E, '-o');
    ylim([-m_E m_E]);
    hold on;
    ylabel("field [V/m]");
    xlabel("z [um]");
    subplot(2,1,2);
    plot(zvec*1e6, sum(P,1), '-o');
    ylim([-m_P m_P]);
    hold on;
    ylabel("polarization density [C/m^2]");
    xlabel("z [um]");
  end
  legend(num2str(dt',"dt = %.0d [s]"));
  title(strcat("ampl = ", num2str(ampl(i))));
  % plot relative error
  figure;
  for j = 1:length(dt)-1
    [E, ~, P, ~] = split_X(X_ref(:,i,j)./p.X_scale, p);
    E_err = abs(E - E_ref)./max(abs(E_ref));
    P_err = abs(sum(P,1) - sum(P_ref,1))./max(abs(sum(P_ref,1)));
    subplot(2,1,1);
    semilogy(zvec*1e6, E_err, '-o');
    ylim([1e-5 1e-1]);
    hold on;
    ylabel("field relative error");
    xlabel("z [um]");
    subplot(2,1,2);
    semilogy(zvec*1e6, P_err, '-o');
    hold on;
    ylim([1e-5 1e-1]);
    ylabel("polarization density relative error");
    xlabel("z [um]");
  end
  legend(num2str(dt(1:end-1)',"dt = %.0d [s]"));
  title(strcat("ampl = ", num2str(ampl(i))));
end
