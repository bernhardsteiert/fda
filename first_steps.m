%% First steps working with live cell imaging data and functional data analysis:
% Warning: The workspace will be cleared!
close all
clear all
clc

fdaMPath = 'fda';
addpath(fdaMPath)

grabdataPath = 'Code + Stage and Outputsignal';
addpath(grabdataPath)

% Stage 17 <--> IGF 100 (from Stage_Treatment_Outputsignal.xlsx)
site = 17;
[timestamp,intensity] = grabdata(site);
c_signal = log10(intensity);

return

%% Reproduce plots from grabdata
close all

plot(timestamp,c_signal,'g','color',[0.7 0.7 0.7])
hold on
plot(timestamp,nanmean(c_signal,2),'color','k','LineWidth',2)

%% Generate spline fits to individual data sets with nbasis basis functions
close all

nbasis = 50;

basis = create_bspline_basis(timestamp(end), nbasis);
smoothed_data = smooth_basis(timestamp,c_signal,basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal,2)))
hold on

plot(smoothed_data)
plot(timestamp,c_signal,'o')

%% Make FPCA - wip
close all

nharm = 2;
c_signal_pcastr = pca_fd(smoothed_data, nharm);

disp(c_signal_pcastr.values(1:4))
plot_pca_fd(c_signal_pcastr, 1, 0)
c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
plot_pca_fd(c_signal_rotpcastr, 1, 0)

