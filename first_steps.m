%% First steps working with live cell imaging data and functional data analysis:
% Warning: The workspace will be cleared!
close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

% Possible stages/sites:
sites = [2, 4, 7, 8, 11, 13, 14, 17, 19, 22, 24, 28, 33, 34, 37, ...
    39, 40, 41, 42, 44, 47, 48, 53, 54, 57, 59, 60, 64, 67, 68];
% Stage 17 <--> IGF 100 (from Stage_Treatment_Outputsignal.xlsx)
% site = 17;
% site = 11; % Unstimulated
site = 2;
if exist(remotepath,'dir')
    [timestamp,intensity] = grabdata(site);
else
    load(['./Workspaces/site_' num2str(site)])
end
c_signal = log10(intensity);

return

%% Reproduce plots from grabdata
close all

plot(timestamp,c_signal,'g','color',[0.7 0.7 0.7])
hold on
plot(timestamp,nanmean(c_signal,2),'color','k','LineWidth',2)

%% Plot every data set with distinct color
close all

% plot_sites = site;
plot_sites = sites;

first_n = 10; % Plot only first_n data-sets

for ip = 1:length(plot_sites)
    if length(plot_sites) > 1
        if exist(remotepath,'dir')
            [timestamp,intensity] = grabdata(plot_sites(ip));
        else
            load(['./Workspaces/site_' num2str(plot_sites(ip))])
        end
        c_signal = log10(intensity);
    end
    
    
    first_n = min(first_n,size(c_signal,2));
    f = figure;
    set(f,'DefaultAxesColorOrder',jet(first_n))

    plot(timestamp,c_signal(:,1:first_n))
    title(['Site ' num2str(plot_sites(ip))])
    
    waitforbuttonpress;
    
    close gcf
end

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

%% Make FPCA with data generated in previous block - wip
close all

nharm = 2;
c_signal_pcastr = pca_fd(smoothed_data, nharm);

disp(c_signal_pcastr.values(1:4))
plot_pca_fd(c_signal_pcastr, 1, 0)
c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Generate smoothing spline fits using as many basis functions as data points
close all

lambda = 1e4; % smoothing parameter
nbasis = length(timestamp);

basis_full = create_bspline_basis(timestamp(end), nbasis);
fdparobj = fdPar(basis_full,2,lambda);
spline_data = smooth_basis(timestamp,c_signal,fdparobj);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal,2)))
hold on

plot(spline_data)
plot(timestamp,c_signal,'o')

%% Make FPCA with data generated in previous block - wip
close all

nharm = 2;
c_signal_pcaspline = pca_fd(spline_data, nharm);

disp(c_signal_pcaspline.values(1:4))
plot_pca_fd(c_signal_pcaspline, 1, 0)
c_signal_rotpcaspline = varmx_pca(c_signal_pcaspline);
plot_pca_fd(c_signal_rotpcaspline, 1, 0)
