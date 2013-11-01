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

% Uncorrelated noise - looks not like measured data (mean too well defined)
close all

noise = .001*randn(8,340);
c_signal = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + noise;

plot(timestamp,c_signal)

% Correlated noise - looks not like measured data (no oscillations)
close all

corr_noise = .0002*randn(8,340);
corr_noise = cumsum(corr_noise,2);
c_signal_corr = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + corr_noise;

plot(timestamp,c_signal_corr)

% Sine modulation
close all

