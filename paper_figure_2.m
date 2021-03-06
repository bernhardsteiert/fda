%% Paper argumentation for fPCA part (figure 2)
% 0. Cut signals in 2 parts dominated by a) Deterministic rise (t = [50 200]) b) Stochastic oscillations (t = [200 650])
% 1. Do fPCA for all ligands at all doses (no inhibitors) for t = [50 200]
%    --> PC2vs3 - IGF one axis - BTC other axis - EGF/HGF/HRG inbetween - EPR/FGF not responding
%    Plot first 4 harmonics --> argument for choosing PCs 2 and 3
% 2. Add EGF inhibitor data and see where it goes (probably: MEKi --> IGF; AKTi --> not responding)
%    OR: Do the same with different inhibitor concentrations
% 3. Subtract mean + first 3 harmonics (plot) from data for t = [200 650]
%    --> Histograms for no dose + highest dose with threshold
%    cell-to-cell heterogeneity <-> time-courses to show oscillations
%    Ligands triggering subpopulation to oscillate; Others don't
% 4. Do same analysis for inhibitor data (EGF and HGF): MEKi induces oscillations; AKTi represses oscillations
%    (Counter-intuitive: EGF-MEKi population mean like IGF, oscillations however not suppressed as for IGF)


close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

% Get properties of sites by calling siteprop(site)
% All ligands + EGF with Inhibitors
sites = [1 2 4:10 17:-1:11 24:30 37:-1:31 41 42 44:50 57:-1:51 64:69];
sites_for_harmonics = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:69];

% sites = [1 4];

times = cell(0);
signals = cell(0);
celltype = [];

for isite = sites
    if exist(remotepath,'dir')
        [times{end+1},intensity] = grabdata(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
        times{end+1} = timestamp;
    end

    log_trafo = 1; % log-transform signal

    if log_trafo
        signals{end+1} = log10(intensity);
    else
        signals{end+1} = intensity;
    end
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);

return

%% Registering of data by calculating mean of pre-stimulation traces
%  Result: Initial behavior does not correlate with late behavior
%% Registering of data by calculating mean of traces at steady-state
% Register inhibitor data from 0 to 35
% t_start  = 0   min
% t_inh    = 40  min
% t_new_ss = 60  min
% t_stim   = 120 min
close all

max_nplot = 50;

t_inh    = 40;
t_new_ss = 60;
t_stim   = 120;
t_ss     = 300;
t_max    = 650;

tmp = find((timestamp - t_inh) < 0);
t_ind_inh = tmp(end);

tmp = find((timestamp - t_new_ss) > 0);
t_ind_new_ss = tmp(1);

tmp = find((timestamp - t_stim) < 0);
t_ind_stim = tmp(end);

tmp = find((timestamp - t_ss) > 0);
t_ind_ss = tmp(1);
[tmp t_ind_max] = min(abs(timestamp - t_max));

c_signal_registered = c_signal;
c_signal_registered2 = c_signal;
c_signal_registered3 = c_signal;

unitype = unique(celltype);

site_lig_ind = [];
site_lig_name = {};
site_lig_dose = [];
site_inh_name = {};
site_inh_dose = [];

n = length(unitype)*4;

rowstocols = 0.6;
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);

ylim = [-.04 .05];
if ~log_trafo
    ylim = 10.^ylim;
end

for iplot = 1:(n./4)
    
    s = siteprop(sites(iplot));
    site_lig_ind = [site_lig_ind s.lig_index];
    site_lig_name{end+1} = s.lig_name;
    site_lig_dose = [site_lig_dose s.lig_dose];
    site_inh_name{end+1} = s.inh_name;
    site_inh_dose = [site_inh_dose s.inh_dose];
    
%     subplot(nrows,ncols,iplot)
%     
%     plot_inds = find(celltype == unitype(iplot));
%     max_plotind = min([length(plot_inds) max_nplot]);
%     
%     plot(timestamp,c_signal(:,plot_inds(1:max_plotind)))
%     hold on
%     plot(timestamp,nanmean(c_signal(:,celltype == unitype(iplot)),2),'k','LineWidth',2)
%     set(gca,'XLim',[0 650])
%     title([site_lig_name{iplot} num2str(site_lig_dose(iplot)) site_inh_name{iplot}])
%     set(gca,'YLim',ylim)
    
    if site_inh_dose(iplot) > 0
        t_ind_start = t_ind_new_ss;
    else
        t_ind_start = 1;
    end
    
    mean_single = nanmean(c_signal_registered(t_ind_start:t_ind_stim,celltype == unitype(iplot)),1);
     
    c_signal_registered(:,celltype == unitype(iplot)) = c_signal(:,celltype == unitype(iplot)) + repmat((nanmean(mean_single) - mean_single),size(c_signal,1),1);
    
%     subplot(nrows,ncols,iplot+2)
% 
%     plot(timestamp,c_signal_registered(:,plot_inds(1:max_plotind)))
%     hold on
%     plot(timestamp,nanmean(c_signal_registered(:,celltype == unitype(iplot)),2),'k','LineWidth',2)
%     set(gca,'XLim',[0 650])
%     title([site_lig_name{iplot} num2str(site_lig_dose(iplot)) site_inh_name{iplot}])
%     set(gca,'YLim',ylim)
%     plot([timestamp(t_ind_start) timestamp(t_ind_start)],ylim,'b--')
%     plot([timestamp(t_ind_stim) timestamp(t_ind_stim)],ylim,'b--')
    
    mean_single2 = nanmean(c_signal_registered2(t_ind_ss:t_ind_max,celltype == unitype(iplot)),1);
     
    c_signal_registered2(:,celltype == unitype(iplot)) = c_signal(:,celltype == unitype(iplot)) + repmat((nanmean(mean_single2) - mean_single2),size(c_signal,1),1);
    
%     subplot(nrows,ncols,iplot+4)
% 
%     plot(timestamp,c_signal_registered2(:,plot_inds(1:max_plotind)))
%     hold on
%     plot(timestamp,nanmean(c_signal_registered2(:,celltype == unitype(iplot)),2),'k','LineWidth',2)
%     set(gca,'XLim',[0 650])
%     title([site_lig_name{iplot} num2str(site_lig_dose(iplot)) site_inh_name{iplot}])
%     set(gca,'YLim',ylim)
%     plot([timestamp(t_ind_ss) timestamp(t_ind_ss)],ylim,'b--')
%     plot([timestamp(t_ind_max) timestamp(t_ind_max)],ylim,'b--')
    
    if site_inh_dose(iplot) > 0
        t_ind_end = t_ind_inh;
    else
        t_ind_end = t_ind_stim;
    end
    
    mean_single3 = nanmean(c_signal_registered3(1:t_ind_end,celltype == unitype(iplot)),1);
     
    c_signal_registered3(:,celltype == unitype(iplot)) = c_signal(:,celltype == unitype(iplot)) + repmat((nanmean(mean_single3) - mean_single3),size(c_signal,1),1);
    
%     subplot(nrows,ncols,iplot+6)
% 
%     plot(timestamp,c_signal_registered3(:,plot_inds(1:max_plotind)))
%     hold on
%     plot(timestamp,nanmean(c_signal_registered3(:,celltype == unitype(iplot)),2),'k','LineWidth',2)
%     set(gca,'XLim',[0 650])
%     title([site_lig_name{iplot} num2str(site_lig_dose(iplot)) site_inh_name{iplot}])
%     set(gca,'YLim',ylim)
%     plot([timestamp(1) timestamp(1)],ylim,'b--')
%     plot([timestamp(t_ind_end) timestamp(t_ind_end)],ylim,'b--')
%     
end

c_signal = c_signal_registered;

%% Register to mean value in time range
time_range = [50 200];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

c_signal = c_signal - repmat(nanmean(c_signal(range_ind,:),1),size(c_signal,1),1);

%% Generate spline fits to data-sets given in sites_for_harmonics
% close all

nbasis = 20;
% time_range = [min(timestamp) max(timestamp)];
time_range = [50 200];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

ind_harm = ismember(celltype,sites_for_harmonics);
ind_fit = ~ind_harm;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_harm),basis);

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

% Define ligands that should be orthogonal
% lig1 = 'IGF';
% lig2 = 'BTC';
% Tricky to do and should be same as set angle by hand:
angle = -10; % rotation angle to right [degree]

Rmat = [1                     0                     0; ...
        0                     cos(2*pi*angle/360)  sin(2*pi*angle/360); ...
        0                     -sin(2*pi*angle/360) cos(2*pi*angle/360)];

flipharm = ones(1,nharm);
flipharm(1:4) = [-1 1 -1 1];

unitypes = unique(celltype(ind_harm));

site_lig_ind = [];
site_lig_name = {};
site_lig_dose = [];
site_inh_name = {};
site_inh_dose = [];

% lig1_ind = [];
% lig2_ind = [];
% lig1_con = [];
% lig2_con = [];

% Exclude Inhibitor data

for isite = sites_for_harmonics   % only plot ligands used for harmonics
    s = siteprop(isite);
    site_lig_ind = [site_lig_ind s.lig_index];
    site_lig_name{end+1} = s.lig_name;
    site_lig_dose = [site_lig_dose s.lig_dose];
    site_inh_name{end+1} = s.inh_name;
    site_inh_dose = [site_inh_dose s.inh_dose];
    
%     if ~isempty(strmatch(s.lig_name,'lig1')) && s.inh_dose == 0
%         lig1_ind = [lig1_ind site_lig_ind];
%         lig1_con = [lig1_con s.lig_dose];
%     end
%     if ~isempty(strmatch(s.lig_name,'lig2')) && s.inh_dose == 0
%         lig2_ind = [lig2_ind site_lig_ind];
%         lig2_con = [lig2_con s.lig_dose];
%     end
end

sites_remain = find(~site_inh_dose);
uni_lig = unique(site_lig_ind(sites_remain));
lig_min = min(site_lig_dose(sites_remain));
lig_max = max(site_lig_dose(sites_remain));
ncolor = 201;
colmap = flipud(jet(ncolor));
color_doses = 10.^linspace(max(log10([lig_min 1])),log10(lig_max),ncolor);

% Rotate to lig1-lig2 basis
% x1 = [x1 repmat(]
% y1 = [y1 flipharm(pcs(1))*c_signal_pcastr.harmscr(celltype(ind_harm) == lig1_ind,pcs(1))];
% y2 = [y2 flipharm(pcs(1))*c_signal_pcastr.harmscr(celltype(ind_harm) == lig2_ind,pcs(2))];

rowstocols = 0.3;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

figure

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

for ilig = 1:length(uni_lig)
    
    subplot(nrows,ncols,ilig)
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
        plot(x_scores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),:),y_scores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),:),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.03) % PC1 vs PC3
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)+.03) % PC1 vs PC2
end

h = get(gca);
subplot(nrows,ncols,ilig+1)
clim = log10([lig_min lig_max]);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([1 2.5 5 10 20 50 100]),'XTickLabel',[0 2.5 5 10 20 50 100])
set(gca,'Visible','off')

%% Figure 2B: Eigenfunctions (new - rotated)
close all
pcs = 1:3;
% Good combination
% angle1 = 30;
% angle2 = -15;

rowstocols = 1;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);

figure

for iplot = 1:nharm
    subplot(nrows,ncols,iplot)
    
    if iplot <= length(pcs)
        tmpplot = sum(repmat(Rmat(iplot,:),size(harm_eval,1),1) .* harm_eval(:,pcs),2);
    else
        tmpplot = harm_eval(:,iplot);
    end

    
    plot(times_fine,tmpplot)
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))]*1.2)
    
    hold on
    plot(time_range,[0 0],'--')
end

%% Figure 2B: Eigenfunctions (old - unrotated)
% close all
angle = 0;

rowstocols = 1;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);

figure

for iplot = 1:nharm
    subplot(nrows,ncols,iplot)
    
    if iplot == pcs(1)
        tmpplot = cos(2*pi*angle/360) * harm_eval(:,iplot) - sin(2*pi*angle/360) * harm_eval(:,pcs(2));
    elseif iplot == pcs(2)
        tmpplot = cos(2*pi*angle/360) * harm_eval(:,iplot) + sin(2*pi*angle/360) * harm_eval(:,pcs(1));
    else
        tmpplot = harm_eval(:,iplot);
    end
    
    plot(times_fine,tmpplot)
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
%     set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))])
    
    hold on
    plot(time_range,[0 0],'--')
end

%% Trying out to rotate eigenfunctions 'by hand' in 3D
close all
pcs = 1:3;
% Good combination
% angle1 = 30;
% angle2 = -15;

angle1 = 25;
angle2 = -20;

R_al = [cos(2*pi*angle1/360)  sin(2*pi*angle1/360)  0; ...
        -sin(2*pi*angle1/360) cos(2*pi*angle1/360)  0; ...
        0                     0                     1];
R_be = [1                     0                     0; ...
        0                     cos(2*pi*angle2/360)  sin(2*pi*angle2/360); ...
        0                     -sin(2*pi*angle2/360) cos(2*pi*angle2/360)];
    
% R_mat = R_be*R_al;
% Reverse order:
R_mat = R_al*R_be;

angle3 = -10;
R_mat = [cos(2*pi*angle3/360) 0                     sin(2*pi*angle3/360); ...
        0                     1                     0; ...
        -sin(2*pi*angle3/360) 0                     cos(2*pi*angle3/360)] * R_mat;

rowstocols = 1;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);

figure

for iplot = 1:nharm
    subplot(nrows,ncols,iplot)
    
    if iplot <= length(pcs)
        tmpplot = sum(repmat(R_mat(iplot,:),size(harm_eval,1),1) .* harm_eval(:,pcs),2);
    else
        tmpplot = harm_eval(:,iplot);
    end

    
    plot(times_fine,tmpplot)
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))]*1.5-.01)
    
    hold on
    plot(time_range,[0 0],'--')
end

%% fPCA with dose-dependend colored ligand level for 'by-hand' rotated harmonics (using R_mat)
% close all
figure

% Define principal components to be plotted
pcs = [1 2];

rowstocols = 0.3;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

for ilig = 1:length(uni_lig)
    
    subplot(nrows,ncols,ilig)
    hold on
    tmpind = find(site_lig_ind(sites_remain) == uni_lig(ilig));
    title(site_lig_name{sites_remain(tmpind(1))})
    
    if ilig == 1
        ylabel(['PC ' num2str(pcs(2))])
    end
    if ilig == length(uni_lig)
        xlabel(['PC ' num2str(pcs(1))])
    end
    
    x_scores = (R_mat(pcs(1),:) * flipped_scores(:,1:3)')';
    y_scores = (R_mat(pcs(2),:) * flipped_scores(:,1:3)')';
    plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);
    
    for isite = tmpind
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_lig_dose(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        plot(x_scores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),:),y_scores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),:),'.','Color',mycolor);
    end
    
%     set(gca,'XLim',[min(cos(2*pi*angle/360)*x_scores-sin(2*pi*angle/360)*y_scores) max(cos(2*pi*angle/360)*x_scores-sin(2*pi*angle/360)*y_scores)]*1.1,'YLim',[min(sin(2*pi*angle/360)*x_scores+cos(2*pi*angle/360)*y_scores) max(sin(2*pi*angle/360)*x_scores+cos(2*pi*angle/360)*y_scores)]*1.1)
    axis equal
end

h = get(gca);
subplot(nrows,ncols,ilig+1)
clim = log10([lig_min lig_max]);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([1 2.5 5 10 20 50 100]),'XTickLabel',[0 2.5 5 10 20 50 100])
set(gca,'Visible','off')

%% Fit additional data with basis from fPCA
% close all

harm_basis = create_fd_basis(c_signal_pcastr.harmfd);
mean_fit = eval_fd(c_signal_pcastr.meanfd,timestamp(range_ind));
smoothed_additional = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_fit)-repmat(mean_fit,1,sum(ind_fit)),harm_basis);

return

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_fit),2)))
hold on


plot(smoothed_additional+c_signal_pcastr.meanfd)
plot(timestamp(range_ind),c_signal(range_ind,ind_fit),'o')

%% Figure 2C: Look what happens with inhibitor data in PCA space
ind_included_panel1 = [4 1 2]; % EGF high dose; + MEKi; + AKTi
ind_included_panel2 = [44 41 42]; % HGF high dose; + MEKi; + AKTi

mycolor = length(ind_included_panel1);

figure

subplot(1,2,1)
plot(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)),flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)),'.','Color',[.7 .7 .7]);
set(gca,'XLim',[min(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1))) max(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2))) max(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)
hold on

xlabel(['PC ' num2str(pcs(1))])
ylabel(['PC ' num2str(pcs(2))])

for isite = ind_included_panel1

    plot(flipharm(pcs(1))*c_signal_pcastr.harmscr(celltype == isite,pcs(1)),flipharm(pcs(2))*c_signal_pcastr.harmscr(celltype == isite,pcs(2)),'.','Color',mycolor);

end