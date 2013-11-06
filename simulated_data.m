%% Try routines on simulated data
close all
% clear all
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

%% Plot: Simulation mean
close all

c_signal_mean = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset;
plot(timestamp,c_signal_mean)

%% Plot: Uncorrelated noise - looks not like measured data (mean too well defined)
close all

noise_level = .001;

noise = noise_level*randn(8,340);
c_signal_sim = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + noise;

plot(timestamp,c_signal_sim)

%% Plot: Correlated noise - looks not like measured data (no oscillations)
close all

noise_level = .0004;

corr_noise = noise_level*randn(8,340);
corr_noise = cumsum(corr_noise,2);
c_signal_corr = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + corr_noise;

plot(timestamp,c_signal_corr)

%% Plot: Sine modulation - looks not like measured data (oscillations too regular)
close all

sine_freq = .01;
mod_amp = 0.1 * mean_amp;
sine_mod = mod_amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) .* sin(2*pi*sine_freq*repmat(timestamp,8,1));
noise_sine = .001*randn(8,340) + sine_mod;
c_signal_sine = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + noise_sine;

plot(timestamp,c_signal_sine)

%% Plot: Use cauchy distributed noise - looks not like measured data (extreme jumps)
noise_level = .0001;

noise_cauchy = noise_level .* randn(8,340) ./ randn(8,340);
c_signal_cauchy = amp./(1+exp(-(repmat(timestamp,8,1)-step_time)/tau_step)) + offset + noise_cauchy;

plot(timestamp,c_signal_cauchy)
set(gca,'YLim',[-.03 .03])

return

%% Plot timecourses from measurements to compare with
figure

plot(timestamp,c_signal)