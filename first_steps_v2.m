%% First steps working with live cell imaging data and functional data analysis:
% Warning: The workspace will be cleared!
% All test code for finding bug in FDA package removed
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
site = 13; % [13 14 17 53 54 57 64 67 68];
if exist(remotepath,'dir')
    [timestamp,intensity] = grabdata(site);
else
    load(['./Workspaces/site_' num2str(site)])
end

log_trafo = 1; % log-transform c_signal

if log_trafo
    c_signal = log10(intensity);
end
range_c_signal = [min(min(c_signal)) max(max(c_signal))];

return

%% Reproduce plots from grabdata
close all

plot(timestamp,c_signal,'g','color',[0.7 0.7 0.7])
hold on
plot(timestamp,nanmean(c_signal,2),'color','k','LineWidth',2)
set(gca,'XLim',[50 650]) % see Pat's Poster

%% Plot every data set with distinct color
close all

plot_sites = site;
% plot_sites = sites;

first_n = 20; % Plot only first_n data-sets

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
    
    set(gca,'XLim',[50 650])
    
    hold on
    plot(timestamp,nanmean(c_signal,2),'color','k','LineWidth',2)
    
    if length(plot_sites) > 1
        waitforbuttonpress;
        close gcf
    end
end

%% Generate spline fits to individual data sets with nbasis basis functions
close all

nbasis = 50;
% time_range = [min(timestamp) max(timestamp)];
time_range = [50 650];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal,2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal(range_ind,:),'o')

%% Check smoothed_data: Reading coefs works, reconstruction of data spline works!
close all

time_range = [50 650];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

tmpcoef = getcoef(smoothed_data);
basis_evaluated = eval_basis(basis,times_fine);

figure
plot(times_fine,basis_evaluated*tmpcoef)
hold on
plot(timestamp(range_ind),c_signal(range_ind,:),'o')

%% Make FPCA with data generated in previous block - wip
close all

% nharm = 6;
nharm = 8;
c_signal_pcastr = pca_fd(smoothed_data, nharm);

% disp(c_signal_pcastr.values(1:4))
plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Plot: %variance explained vs. #basis functions
close all

thres_var = 0.9;

cumprobs = cumsum([0;c_signal_pcastr.varprop]);
[tmp thres_ind] = min(abs(cumprobs - thres_var));

plot(0:length(c_signal_pcastr.varprop),cumprobs)
hold on
plot(0:thres_ind-1,ones(1,thres_ind)*thres_var,'--')
plot([thres_ind-1 thres_ind-1],[0 thres_var],'--')

xlabel('fPCA basis functions')
ylabel('cumulative variance explained')

fprintf('To explain %s variance, use %i fPCA basis functions.\n\n',num2str(thres_var,3),thres_ind-1);

%% Alternative: Use pareto function of Matlab
close all

thres_var = 0.9;

percent_explained = 100*c_signal_pcastr.varprop;

pareto(percent_explained)
% hold on
% plot(0:thres_ind-1,ones(1,thres_ind)*100*thres_var,'--')
% plot([thres_ind-1 thres_ind-1],[0 100*thres_var],'--')

xlabel('Principal Component')
ylabel('Variance Explained (%)')

%% Plot: Harmonic scores PC1 vs. PC2
close all

plot(c_signal_pcastr.harmscr(:,1),c_signal_pcastr.harmscr(:,2),'.')

%% Plot harmonics from harmfd object
close all

nperplot = 4;
nsubplots = ceil(size(c_signal,2)./nperplot);
rowstocols = 0.5;
time_range = [50 650];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);
harm_fine = eval_fd(c_signal_pcastr.harmfd,times_fine);

nrows = ceil(nsubplots^rowstocols);
ncols = ceil(nsubplots / nrows);

data_fpca_repr = c_signal_pcastr.harmscr*harm_fine';
mean_fine = eval_fd(c_signal_pcastr.meanfd,times_fine);

for iplot = 1:nsubplots
    subplot(nrows,ncols,iplot)
    
    inds = ((iplot-1)*nperplot+1):min([size(c_signal,2) iplot*nperplot]);
    plot(times_fine,repmat(mean_fine,1,length(inds))+data_fpca_repr(inds,:)')
    hold on
    plot(timestamp(range_ind),c_signal(range_ind,inds),'o')

end
    
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

%% Plot data with mean subtracted
close all

f = figure;
set(f,'DefaultAxesColorOrder',hsv(first_n))

plot(timestamp,c_signal-repmat(nanmean(c_signal,2),1,size(c_signal,2)))

set(gca,'XLim',[50 650])
