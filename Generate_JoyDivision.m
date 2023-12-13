% %% HIGH FIELD SINUSOIDAL
% clear
% close all
% 
% [p, newton_opts, trap_opts] = demo_initialize(40e-6, 1001);
% 
% p.ampl_J = 4e7/p.dz*p.omega_J;  % source amplitude
% p.dt = 1e-16;                   % minimum timestep
% p.tf = 110e-15;                  % final simulation time
% p.source_type = "sinusoid";       % source type
% 
% % initial conditions
% E0 = zeros(p.N,1);
% dtE0 = zeros(p.N, 1);
% P0 = zeros(size(p.Lorentz,1), p.N);
% dtP0 = zeros(size(p.Lorentz,1), p.N);
% X0 = generate_X(E0, dtE0, P0, dtP0, p);
% trap_opts.visualize_dt = Inf;
% 
% tic
% % compute preconditioner
% Jf_sin_high = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% % trapezoidal method
% [t_sin_high, X_sin_high] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf_sin_high, trap_opts, newton_opts);
% toc
% 
% %% FFT SINUSOIDAL
% movie_fft(X_sin_high, t_sin_high, p, false)
% 
% %% COLORMAP
% space_colormap(X_sin_high, t_sin_high, p)

%% LOW FIELD
clear
close all

[p, newton_opts, trap_opts] = demo_initialize(7e-6, 501);

p.ampl_J = 1e7/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-16;                   % minimum timestep
p.tf = 20e-15;                  % final simulation time
p.source_type = "ricker";       % source type

% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = 2e-16;
trap_opts.visualize_dt = Inf;

% compute preconditioner
Jf = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
% trapezoidal method
[t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf, trap_opts, newton_opts);

%%
[E, ~, ~, ~] = split_X(X(:, end), p);
plot(E)
center = 369;
width = 84;
points = (center-width):(center+width);


z_temp = 0:p.dz:(p.dz*size(E, 1)-p.dz);
z = z_temp - mean(z_temp);

plot(z(points), E(points))
xlim([z(points(1)), z(points(end))])


%%
clear
close all

amplitudes = logspace(6, 8, 20);

[p, newton_opts, trap_opts] = demo_initialize(7e-6, 501);
p.ampl_J = 1e7/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-16;                   % minimum timestep
p.tf = 20e-15;                  % final simulation time
p.source_type = "ricker";       % source type
% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = Inf;
center = 369;
width = 84;
points = (center-width):(center+width);
z_temp = 0:p.dz:(p.dz*501-p.dz);
z = z_temp - mean(z_temp);


for amp_i = amplitudes
    p.ampl_J = amp_i/p.dz*p.omega_J;  % source amplitude
    Jf = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
    [t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf, trap_opts, newton_opts);

    [E, ~, ~, ~] = split_X(X(:, end), p);
    E = E/norm(E);
    plot(z(points), E(points))
    xlim([z(points(1)), z(points(end))])
    pause(0.1)
end

%%

clear
close all

tic
amplitudes = logspace(6, 8.1, 50);

[p, newton_opts, trap_opts] = demo_initialize(7e-6, 501);
p.ampl_J = 1e7/p.dz*p.omega_J;  % source amplitude
p.dt = 1e-16;                   % minimum timestep
p.tf = 20e-15;                  % final simulation time
p.source_type = "ricker";       % source type
% initial conditions
E0 = zeros(p.N,1);
dtE0 = zeros(p.N, 1);
P0 = zeros(size(p.Lorentz,1), p.N);
dtP0 = zeros(size(p.Lorentz,1), p.N);
X0 = generate_X(E0, dtE0, P0, dtP0, p);
trap_opts.visualize_dt = Inf;
center = 355;
% width = 64;
width = 105;
points = (center-width):(center+width);
z_temp = 0:p.dz:(p.dz*501-p.dz);
z = z_temp - mean(z_temp);
E_vect = zeros(length(points), length(amplitudes));

% Define the Brewer colormap
colormap_brewer = brewermap(length(amplitudes), 'RdBu');
spacing = 0.025;
k = 1;

figure;

for amp_i = 1:length(amplitudes)
    p.ampl_J = amplitudes(amp_i)/p.dz*p.omega_J;  % source amplitude
    Jf = FastJacobian(@(x) eval_f(x, p, 0), X0, p, 1e-6, p.N);
    [t,X] = trapezoid(@eval_f, p, @eval_u, X0, p.tf, p.dt, Jf, trap_opts, newton_opts);

    [E, ~, ~, ~] = split_X(X(:, end), p);
    E = E/norm(E);
    
    % Calculate the vertical position for plotting
    vertical_position = (amp_i - 1) * spacing;

    % Plot the function with the Brewer colormap
    % plot(z(points), E(points) + vertical_position, 'Color', colormap_brewer(amp_i, :), 'linewidth', 2);
    plot(z(points), E(points) + vertical_position, 'k', 'linewidth', 2);
    E_vect(:, k) = E(points);
    k = k+1;
    hold on;
end

z_vect = z(points);

% Customize the plot
xlabel('Position along z-axis');
ylabel('Function Value');
title('Functions with Brewer Colormap');
colorbar('Ticks', linspace(0, 1, 40), 'TickLabels', num2str(amplitudes', '%.2e'));
colormap(colormap_brewer);
% caxis([min(amplitudes), max(amplitudes)]);

hold off;
toc

%%

tic
scala = 10;
% Define the Brewer colormap
colormap_brewer = brewermap(length(amplitudes)+scala, 'Blues');
spacing = 0.025;
margin = 0.3e-6;

fft_fig = figure(1);
set(fft_fig, 'resize', 'off', 'Position', [100 100 800 600]);
z_vect_0 = z_vect-min(z_vect);

for i = 1:length(E_vect(1, :))
    
    vertical_position = (i - 1) * spacing;
    % plot(z_vect, E_vect(:, i) + vertical_position, 'k', 'linewidth', 2);
    plot(z_vect_0, E_vect(:, i) + vertical_position, 'Color', colormap_brewer(i+scala, :), 'linewidth', 2);
    hold on;
end

% Customize the plot
xlabel('Position (m)');
ylabel('Electric field (a. u)');
% title('Functions with Brewer Colormap');
xlim([min(z_vect_0)-margin, max(z_vect_0)+margin])
ax = gca;
ax.FontSize = 18;
hold off;
toc

%%

tic
scala = 10;
% Define the Brewer colormap
colormap_brewer = brewermap(length(amplitudes)+scala, 'Blues');
spacing = 0.025;
margin = 0.3e-6;

fft_fig = figure(1);
set(fft_fig, 'resize', 'off', 'Position', [100 100 400 600]);

% Set the figure background color to black
set(gcf, 'Color', 'k');

z_vect_0 = z_vect-min(z_vect);

for i = 1:length(E_vect(1, :))
    vertical_position = (i - 1) * spacing;
    
    % Plot the curves with light grey color on black background
    % plot(z_vect_0, E_vect(:, i) + vertical_position, 'Color', [0.8 0.8 0.8], 'linewidth', 2);
    plot(z_vect_0, E_vect(:, i) + vertical_position, 'Color', [1 1 1], 'linewidth', 2);
    hold on;
end

% Customize the plot
xlabel('Position (m)', 'Color', 'w'); % Set xlabel color to white
ylabel('Electric field (a. u)', 'Color', 'w'); % Set ylabel color to white
xlim([min(z_vect_0)-margin, max(z_vect_0)+margin])
ax = gca;
ax.FontSize = 18;
ax.Color = 'k'; % Set axis background color to black

% Set grid color to white
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca,'GridColor','w');

set(gcf, 'InvertHardcopy', 'off', 'Color', 'k');
% print('joydivision.png', '-dpng', '-r900'); % Change the filename as needed

hold off;
toc
