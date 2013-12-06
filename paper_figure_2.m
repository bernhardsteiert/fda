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

return

plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Figure 2A: fPCA with dose-dependend colored ligand level
% close all

% Define principal components to be plotted
pcs = [2 3];

flipharm = ones(1,nharm);
flipharm(1:4) = [1 1 -1 1];

unitypes = unique(celltype(ind_harm));

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

sites_remain = find(~site_inh_dose);
uni_lig = unique(site_lig_ind(sites_remain));
lig_min = min(site_lig_dose(sites_remain));
lig_max = max(site_lig_dose(sites_remain));
ncolor = 201;
colmap = flipud(jet(ncolor));
color_doses = 10.^linspace(max(log10([lig_min 1])),log10(lig_max),ncolor);

rowstocols = 0.3;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

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
    
    plot(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)),flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)),'.','Color',[.7 .7 .7]);
    
    for isite = tmpind
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_lig_dose(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        plot(flipharm(pcs(1))*c_signal_pcastr.harmscr(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),pcs(1)),flipharm(pcs(2))*c_signal_pcastr.harmscr(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),pcs(2)),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1))) max(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2))) max(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)
    
end

h = get(gca);
subplot(nrows,ncols,ilig+1)
clim = log10([lig_min lig_max]);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([1 2.5 5 10 20 50 100]),'XTickLabel',[0 2.5 5 10 20 50 100])
set(gca,'Visible','off')

%% Figure 2B: Eigenfunctions
% close all

rowstocols = 1;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);

for iplot = 1:nharm
    subplot(nrows,ncols,iplot)
    
    plot(times_fine,harm_eval(:,iplot))
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))])
    
    hold on
    plot(time_range,[0 0],'--')
end

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