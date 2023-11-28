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

if ~isfile("reference_accuracy_ampl_5e7.mat");
  opts.visualize = Inf;
  opts.save_intermediate = false;
  dt = logspace(-19,-21,3);
  X_ref = zeros(size(X0,1),length(dt));
  ref_confidence = zeros(length(dt));
  disp("computing reference solutions X_ref");
  for i = 1:length(dt)
    disp(strcat("computing reference for timestep dt = ", num2str(dt(i)*1e18), " [as]"));
    tic;
    [t,X] = forward_euler(@eval_f, p, @eval_u, X0, p.tf, dt(i), opts);
    toc;
    X_ref(:,i) = X;
    if i == 1
      ref_confidence(i) = Inf;
    else
      ref_confidence(i) = max(abs(X - X_ref(:,i-1)));
    end
    save("reference_accuracy_ampl_5e7.mat","X_ref","ref_confidence","dt");
  end
else
  load("reference_accuracy_ampl_5e7.mat");
  % X_ref, "ref_confidence", "dt"
end

figure;
zvec = linspace(0, (p.N-1)*p.dz, p.N);
for i = 1:length(dt)
  [E, ~, P, ~] = split_X(X_ref(:,i)./p.X_scale, p);
  subplot(2,1,1);
  plot(zvec*1e6, E, '-o');
  hold on;
  ylabel("field [V/m]");
  xlabel("z [um]");
  subplot(2,1,2);
  plot(zvec*1e6, sum(P,1), '-o');
  hold on;
  ylabel("polarization density [C/m^2]");
  xlabel("z [um]");
end
legend("dt = 10^{-19} [s]", "dt = 10^{-20} [s]", "dt = 10^{-21} [s]");
% plot relative error
[E_ref, ~, P_ref, ~] = split_X(X_ref(:,length(dt))./p.X_scale, p);
for i = 1:length(dt)-1
  [E, ~, P, ~] = split_X(X_ref(:,i)./p.X_scale, p);
  E_err = abs(E - E_ref)./max(abs(E_ref));
  P_err = abs(sum(P,1) - sum(P_ref,1))./max(abs(sum(P_ref,1)));
  subplot(2,1,1);
  plot(zvec*1e6, E_err, '-o');
  hold on;
  ylabel("field relative error");
  xlabel("z [um]");
  subplot(2,1,2);
  plot(zvec*1e6, P_err, '-o');
  hold on;
  ylabel("polarization density relative error");
  xlabel("z [um]");
end
legend("dt = 10^{-19} [s]", "dt = 10^{-20} [s]", "dt = 10^{-21} [s]");
