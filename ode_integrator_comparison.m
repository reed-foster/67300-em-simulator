clear all
close all

p = simulation_params();

% override timestep
% p.dt = 1e-20;
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
% t_unstable > 1e-19 for forward_euler at all field strengths

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

% evaluate forward euler
for i=1:length(ampl)
  figure;
  zvec = linspace(0, (p.N-1)*p.dz, p.N);
  [E_ref, ~, P_ref, ~] = split_X(X_ref(:,i,end)./p.X_scale, p);
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

% now compare Forward Euler with Trapezoidal + NewtonGCR
% sweep over newton err_rel and err_gcr
% also sweep over timestep and excitation amplitude
% run each sim 50x and pick the median of the top 10 times
dt_euler_unstable = 1e-19;
N_trials = 20;
dt = logspace(-15, -19, 13);
ampl = logspace(-4, 8, 4);
t_euler = zeros(N_trials, length(ampl), length(dt));
t_trap = zeros(N_trials, length(ampl), length(dt));
t_trap_precond = zeros(N_trials, length(ampl), length(dt));
X_euler = zeros(size(X0,1), length(ampl), length(dt));
X_trap = zeros(size(X0,1), length(ampl), length(dt));
X_trap_precond = zeros(size(X0,1), length(ampl), length(dt));

% setup for Trapezoidal + NewtonGCR
newton_opts.err_f = Inf;
newton_opts.err_dx = Inf;
newton_opts.err_rel = 1e-2;
newton_opts.max_iter = 20;
newton_opts.matrix_free = true; % use matrix-free GCR solver for dx = J\(-f)
newton_opts.err_gcr = 1e-2; % relative error for GCR residual inside Newton
newton_opts.eps_fd = 1e-7; % relative perturbation for Jacobian
trap_opts.save_intermediate = false; % only save final state
trap_opts.visualize_dt = Inf; % don't visualize
trap_opts.print_debug = false; % don't print number of iterations
trap_opts.adaptive_timestep = false;
trap_opts.linear_only = false;
euler_opts.save_intermediate = false;
euler_opts.visualize_dt = Inf;

% zero-field Jacobian
[Jf0, ~] = Linearize_f(@eval_f, X0, p, zeros(p.N,1), 1e-6);

if ~isfile("ode_comparison.mat")
  for trial=1:N_trials
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    disp(num2str(trial,"running trial %d"));
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    for ampl_i=1:length(ampl)
      p.ampl_J = ampl(ampl_i)/p.dz*p.omega_J;
      disp(num2str(p.ampl_J,"setting amplitude to %.2d [A/m/s]"));
      for dt_i=1:length(dt)
        p.dt = dt(dt_i);
        disp(num2str(p.dt,"setting timestep to %.2d [s]"));
        % do euler measurement
        disp("doing Forward Euler measurement");
        tic;
        [t,X] = forward_euler(@eval_f, p, @eval_u, X0, p.tf, p.dt, euler_opts);
        t_euler(trial,ampl_i,dt_i) = toc;
        X_euler(:,ampl_i,dt_i) = X;
        disp(num2str(t_euler(trial,ampl_i,dt_i), "... took %d [s]"));
        save("ode_comparison.mat","N_trials","ampl","dt","t_euler","t_trap","t_trap_precond","X_euler","X_trap","X_trap_precond","newton_opts","trap_opts","p","Jf0");

        if p.dt > 5e-19
          % do trap measurement
          disp("doing Trap + Newton GCR measurement (no preconditioner)");
          newton_opts.preconditioner = false;
          tic;
          [t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf0, trap_opts, newton_opts);
          t_trap(trial,ampl_i,dt_i) = toc;
          X_trap(:,ampl_i,dt_i) = X;
          disp(num2str(t_trap(trial,ampl_i,dt_i), "... took %d [s]"));

          disp("doing Trap + Newton GCR measurement (with preconditioner)");
          newton_opts.preconditioner = true;
          tic;
          [t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf0, trap_opts, newton_opts);
          t_trap_precond(trial,ampl_i,dt_i) = toc;
          X_trap_precond(:,ampl_i,dt_i) = X;
          disp(num2str(t_trap_precond(trial,ampl_i,dt_i), "... took %d [s]"));
          save("ode_comparison.mat","N_trials","ampl","dt","t_euler","t_trap","t_trap_precond","X_euler","X_trap","X_trap_precond","newton_opts","trap_opts","p","Jf0");
        end
      end
    end
  end
else
  load("ode_comparison.mat");
  % N_trials, ampl, dt, err_rel, err_gcr, t_euler, t_trap, X_euler, X_trap, newton_opts, trap_opts, p, Jf0
end

% compare
% do plot of error vs runtime with separate series for each amplitude
t_euler_avg = [];
t_trap_avg = [];
t_trap_precond_avg = [];
err_euler = [];
err_trap = [];
err_trap_precond = [];
N_avg = 10;
f1 = figure; % low amplitude plot
f2 = figure; % high amplitude plot
C = colororder;
field_ampl = [];
for ampl_i = [1 length(ampl)]
  [E_ref, ~, P_ref, ~] = split_X(X_ref(:,ampl_i,end)./p.X_scale, p);
  field_ampl(end+1) = max(abs(E_ref));
  for dt_i=1:length(dt)
    % get error and avg runtime of Forward Euler
    t_euler_avg(ampl_i,end+1) = median(mink(t_euler(:,ampl_i,dt_i),N_avg));
    [E, ~, P, ~] = split_X(X_euler(:,ampl_i,dt_i)./p.X_scale, p);
    E_err = max(abs(E - E_ref)./max(abs(E_ref)));
    P_err = max(abs(sum(P,1) - sum(P_ref,1))./max(abs(sum(P_ref,1))));
    err_euler(ampl_i,end+1) = max(E_err, P_err);

    % get error and avg runtime of Trapezoidal (without preconditioner)
    runtime = median(mink(t_trap(:,ampl_i,dt_i),N_avg));
    if runtime > 0
      [E, ~, P, ~] = split_X(X_trap(:,ampl_i,dt_i)./p.X_scale, p);
      E_err = max(abs(E - E_ref)./max(abs(E_ref)));
      P_err = max(abs(sum(P,1) - sum(P_ref,1))./max(abs(sum(P_ref,1))));
      if (max(E_err, P_err) <= 1)
        err_trap(ampl_i,end+1) = max(E_err, P_err);
        t_trap_avg(ampl_i,end+1) = runtime;
      end
    end

    % get error and avg runtime of Trapezoidal (with preconditioner)
    runtime = median(mink(t_trap_precond(:,ampl_i,dt_i),N_avg));
    if runtime > 0
      [E, ~, P, ~] = split_X(X_trap_precond(:,ampl_i,dt_i)./p.X_scale, p);
      E_err = max(abs(E - E_ref)./max(abs(E_ref)));
      P_err = max(abs(sum(P,1) - sum(P_ref,1))./max(abs(sum(P_ref,1))));
      if (max(E_err, P_err) <= 1)
        err_trap_precond(ampl_i,end+1) = max(E_err, P_err);
        t_trap_precond_avg(ampl_i,end+1) = runtime;
      end
    end
  end
  if ampl_i == 1
    figure(f1);
  else
    figure(f2);
  end
  h = loglog(t_euler_avg(ampl_i,:), err_euler(ampl_i,:), '^', MarkerSize=10);
  set(h, 'MarkerFaceColor', get(h, 'Color'));
  hold on;
  h = loglog(t_trap_avg(ampl_i,:), err_trap(ampl_i,:), 'o', MarkerSize=10);
  set(h, 'MarkerFaceColor', get(h, 'Color'));
  h = loglog(t_trap_precond_avg(ampl_i,:), err_trap_precond(ampl_i,:), 's', MarkerSize=10);
  set(h, 'MarkerFaceColor', get(h, 'Color'));
  legend('Forward Euler', 'Trapezoidal/Newton-GCR', 'Trapezoidal/Newton-GCR + Preconditioner');
  ylim([1e-3 1e1]);
  xlim([0.5 50]);
end

figure(f1);
ax = gca;
ax.FontSize = 16;
title(num2str(field_ampl(1), "E-field amplitude = %.0e [V/m]"));
xlabel("average runtime [s]",FontSize=16);
ylabel("relative error of output (L_\infty)",FontSize=16);

figure(f2);
ax = gca;
ax.FontSize = 16;
title(num2str(field_ampl(end), "E-field amplitude = %.0e [V/m]"));
xlabel("average runtime [s]",FontSize=16);
ylabel("relative error of output (L_\infty)",FontSize=16);

