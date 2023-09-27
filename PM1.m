clear all
close all
clc

% nodal variables
N = 500; % discretization
E = zeros(N); % electric field 
E_dt = zeros(N); % first derivative of E
H = zeros(N); % magnetic field 
P = zeros(N); % polarization vector
P_dt = zeros(N); % first derivative of P

%parameters 
eps_0 =  8.85e-12; % F/m
mu_0 = 1.26e-6; % N/A^2
omega_chi_1 = 2*pi*2.9e15; % Hz
chi_1 = 4;
chi_2 = 41.7e-12; % m/V %Laboratory for Nanoscale Optics, John A. Paulson School of Engineering and Applied Sciences, Harvard University
delta_chi_1 = 2*pi*1.2e13; % Hz
delta_x = 50e-9; % m

% source 
J = zeros(N);


