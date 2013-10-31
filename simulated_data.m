%% Try routines on simulated data
close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

%% Generate simulations by hand
% step-like behavior for site = 17 (IGF high-dose)
logyrange = [-.015 .015];
timestamp = linspace(0,1.6981e+03,340);

% step approximately at 127.9 min
step_time = 127.9;
[m step_ind] = min(abs(timestamp - step_time));
tau_step = 8;
mean_amp = range(logyrange);
offset = -.015;

amp_randn = randn(8,1)*ones(1,340);
amp = mean_amp + .15*mean_amp*amp_randn;

off_randn = randn(8,1)*ones(1,340);
offset = offset + .2*offset*off_randn;

% Uncorrelated noise, better use correlated noise or modulate with sine
noise = .001*randn(8,340);
c_signal = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + noise;

plot(timestamp,c_signal)