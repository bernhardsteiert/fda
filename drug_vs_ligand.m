%% Make complete fPCA for data without drugs and drug / EGF concentration dependend

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

log_trafo = 1; % log-transform signal
time_range = [50 200];

% LIGAND DATA

% Get properties of sites by calling siteprop(site)
sites_for_harmonics = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:69];

% sites = [1 4];

times = cell(0);
signals = cell(0);
celltype = [];

for isite = sites_for_harmonics
    if exist(remotepath,'dir')
        [times{end+1},intensity] = grabdata(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
        times{end+1} = timestamp;
    end

    if log_trafo
        signals{end+1} = log10(intensity);
    else
        signals{end+1} = intensity;
    end
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

% DRUG DATA

dataPath = '2D_dose_response_drugsVSEGF_130903';

% site 20 (row 2 - col 1) is missing ...
sites_drug = [1:19 21:70];
lig_name = 'EGF';

times_drug = cell(0);
signals_drug = cell(0);
celltype_drug = [];

for isite = sites_drug
    if exist(remotepath,'dir')
        [times_drug{end+1},intensity] = grabdata(isite,dataPath);
    else
        load(['./Workspaces/site_' num2str(isite) '_only_' lig_name])
        times_drug{end+1} = timestamp;
    end

    if log_trafo
        signals_drug{end+1} = log10(intensity);
    else
        signals_drug{end+1} = intensity;
    end
    
    celltype_drug = [celltype_drug ones(1,size(intensity,2))*isite];
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);

timestamp_drug = times_drug{1}; % same time sampling for all data sets
timestamp_drug = timestamp_drug + 10; % Timing here is 10 min shifted to timing in first data-sets
c_signal_drug = cell2mat(signals_drug);

return

%% Plot raw data of defined ligands
close all

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

figure

posFig = get(gcf,'Position');
% posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./10);

site_lig_ind = [];
site_lig_name = {};
site_lig_dose = [];
site_inh_name = {};
site_inh_dose = [];

% Exclude Inhibitor data
for isite = sites_for_harmonics   % only plot ligands used for harmonics
    s = siteprop(isite);
    site_lig_ind = [site_lig_ind s.lig_index];
    site_lig_name{end+1} = s.lig_name;
    site_lig_dose = [site_lig_dose s.lig_dose];
    site_inh_name{end+1} = s.inh_name;
    site_inh_dose = [site_inh_dose s.inh_dose];
end

nrows = 7;
ncols = 7;

for iplot = 1:length(sites_for_harmonics)
    subplot(nrows,ncols,iplot)
    
    plot(repmat(timestamp(range_ind),1,sum(celltype == sites_for_harmonics(iplot))),c_signal(range_ind,celltype == sites_for_harmonics(iplot)),'g','color',[0.7 0.7 0.7])
    hold on
    plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == sites_for_harmonics(iplot)),2),'color','k','LineWidth',2)
    title([site_lig_name{iplot} num2str(site_lig_dose(iplot))])
    
    ylim = [-1 1]*.04;
    if ~log_trafo
        ylim = 10.^ylim;
    end
    set(gca,'XLim',time_range,'YLim',ylim)

end


%% Reproduce plot from Bernhard Kraemer
% close all

[tmp range_ind_min] = min(abs(timestamp_drug - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp_drug - time_range(2)));
range_ind_drug = range_ind_min:range_ind_max;

site_lig_dose_drug = [];
site_inh_ind_drug = [];
site_inh_name_drug = {};
site_inh_dose_drug = [];

for isite = sites_drug
    s = siteprop_drug(isite);
    site_lig_dose_drug = [site_lig_dose_drug s.lig_dose];
    site_inh_ind_drug = [site_inh_ind_drug s.inh_ind];
    site_inh_name_drug{end+1} = s.inh_name;
    site_inh_dose_drug = [site_inh_dose_drug s.inh_dose];
end

figure

posFig = get(gcf,'Position');
% posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./8);


nrows = 7;
ncols = 10;

for iplot = sites_drug
    isite = iplot;
    % Switch every second row
    inh_index = ceil(iplot/10);
    dose_index = mod(iplot-1,10)+1;
    if ~mod(inh_index,2)
        dose_index = 11-dose_index;
    end
    isite = 10*(inh_index-1)+dose_index;
    
    subplot(nrows,ncols,isite)
    hold on
    
    plot(repmat(timestamp_drug(range_ind_drug),1,sum(celltype_drug == iplot)),c_signal_drug(range_ind_drug,celltype_drug == iplot),'g','color',[0.7 0.7 0.7])
    plot(timestamp_drug(range_ind_drug),nanmean(c_signal_drug(range_ind_drug,celltype_drug == iplot),2),'color','k','LineWidth',2)
    prop_ind = (iplot == sites_drug);
    title([lig_name ' ' num2str(site_lig_dose_drug(prop_ind)) site_inh_name_drug{prop_ind} num2str(site_inh_dose_drug(prop_ind))])
    
    ylim = [-1 1]*.02;
    if ~log_trafo
        ylim = 10.^ylim;
    end
    set(gca,'XLim',time_range,'YLim',ylim)

end


%% Generate spline fits to data-sets given in sites_for_harmonics and sites_drug
% close all

nbasis = 20;

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

[tmp range_ind_min] = min(abs(timestamp_drug - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp_drug - time_range(2)));
range_ind_drug = range_ind_min:range_ind_max;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data_lig = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),basis);

smoothed_data_drug = smooth_basis(timestamp_drug(range_ind_drug),c_signal_drug(range_ind_drug,:),basis);

smoothed_data = horzcat(smoothed_data_lig,smoothed_data_drug);

return

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_harm),2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal(range_ind,ind_harm),'o')

%% Make FPCA with data generated in previous block
% close all

nharm = 4;
% c_signal_pcastr = pca_fd(smoothed_data, nharm);
c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!
% c_signal_pcastr = varmx_pca(c_signal_pcastr);

return

plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Figure 2A: fPCA with dose-dependend colored ligand level
close all

% Define principal components to be plotted
pcs = [2 3];

% angle = 0;
angle = -10; % rotation angle to right [degree]

Rmat = [1                     0                     0; ...
        0                     cos(2*pi*angle/360)  sin(2*pi*angle/360); ...
        0                     -sin(2*pi*angle/360) cos(2*pi*angle/360)];

flipharm = ones(1,nharm);
flipharm(1:4) = [-1 1 -1 1];

unitypes = unique(celltype);

sites_remain = 1:length(sites_for_harmonics);
uni_lig = unique(site_lig_ind(sites_remain));
lig_min = min(site_lig_dose(sites_remain));
lig_max = max(site_lig_dose(sites_remain));
ncolor = 201;
colmap = flipud(jet(ncolor));
color_doses = 10.^linspace(max(log10([lig_min 1])),log10(lig_max),ncolor);

rowstocols = 0.3;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

figure

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

sps = [];

for ilig = 1:length(uni_lig)
    
    tmp = subplot(nrows,ncols,ilig);
    sps = [sps tmp];
    hold on
    tmpind = find(site_lig_ind(sites_remain) == uni_lig(ilig));
    title(site_lig_name{sites_remain(tmpind(1))})
    
    if ilig == 1
        ylabel(['PC ' num2str(pcs(2))])
    end
    if ilig == length(uni_lig)
        xlabel(['PC ' num2str(pcs(1))])
    end
    
    x_scores = (Rmat(pcs(1),:) * flipped_scores(:,1:3)')';
    y_scores = (Rmat(pcs(2),:) * flipped_scores(:,1:3)')';
    plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);
    
    for isite = tmpind
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_lig_dose(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        plot(x_scores(celltype == sites_for_harmonics(sites_remain(isite)),:),y_scores(celltype == sites_for_harmonics(sites_remain(isite)),:),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.03) % PC1 vs PC3
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)+.03) % PC1 vs PC2
    set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.1) % PC2 vs PC3
end

set(sps(1),'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.1) % PC2 vs PC3

h = get(gca);
subplot(nrows,ncols,ilig+1)
clim = log10([lig_min lig_max]);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([1 2.5 5 10 20 50 100]),'XTickLabel',[0 2.5 5 10 20 50 100])
set(gca,'Visible','off')


% DRUG CASE

figure

sites_remain = 1:length(sites_drug);
uni_lig = unique(site_lig_dose_drug);
uni_drug = unique(site_inh_ind_drug(sites_remain));
drug_min = min(site_inh_dose_drug(sites_remain));
drug_max = max(site_inh_dose_drug(sites_remain));
ncolor = 201;
colmap = flipud(jet(ncolor));
color_doses = 10.^linspace(max(log10([drug_min .1])),log10(drug_max),ncolor);

rowstocols = 0.3;
nrows = length(uni_drug);
ncols = length(uni_lig)+1;

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/3;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./10);
sps = [];

for idrug = 1:length(uni_drug)
    for ilig = 1:length(uni_lig)
    
    tmp = subplot(nrows,ncols,(idrug-1)*ncols+ilig);
    sps = [sps tmp];
    hold on
    tmpind = find(site_inh_ind_drug(sites_remain) == uni_drug(idrug));
    tmpind2 = find(site_lig_dose_drug(sites_remain) == uni_lig(ilig));
    title([lig_name ' '  num2str(site_lig_dose_drug(tmpind2(1))) site_inh_name_drug{sites_remain(tmpind(1))}])
    
    if idrug*ilig == 1
        ylabel(['PC ' num2str(pcs(2))])
    end
    if idrug == length(uni_drug) && ilig == 1
        xlabel(['PC ' num2str(pcs(1))])
    end
    
    x_scores = (Rmat(pcs(1),:) * flipped_scores(:,1:3)')';
    y_scores = (Rmat(pcs(2),:) * flipped_scores(:,1:3)')';
    plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);
    
    for isite = intersect(tmpind,tmpind2)
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_inh_dose_drug(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        drug_ind = size(c_signal,2) + find(celltype_drug == sites_drug(sites_remain(isite)));
        plot(x_scores(drug_ind,:),y_scores(drug_ind,:),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.03) % PC1 vs PC3
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)+.03) % PC1 vs PC2
    set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.1) % PC2 vs PC3
    end
end

set(sps(1),'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.1) % PC2 vs PC3

h = get(gca);
subplot(nrows,ncols,[((idrug-1)*ncols+ilig+1)./2 (idrug-1)*ncols+ilig+1])
clim = log10([drug_min drug_max]);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(flipud(jet(ncolor)))
colorbar('YTick',log10(linspace(1,10,7)),'YTickLabel',[0 0.25 0.5 1 2.5 5 10])
set(gca,'Visible','off')
