%% Make analysis of pulsatory behavior for data without drugs and drug / EGF concentration dependend
% BEWARE: Signals still not scaled between experiments
% If shift of pulsatory behavior with rising inhibitor is observed, it should still be OK ...

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

log_trafo = 1; % log-transform signal
time_range = [200 515];

% Gates for drug signals
% lower_gate_all = 0.97;
% lower_gate_late = 0.98;
lower_gate_all = 0.93;
lower_gate_late = 0.95;

if log_trafo
    lower_gate_all = log10(lower_gate_all);
    lower_gate_late = log10(lower_gate_late);
end

% LIGAND DATA

% Get properties of sites by calling siteprop(site)
sites_for_harmonics = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:69];

% sites = [1 4];

times = cell(0);
signals = cell(0);
signals_raw = cell(0);
celltype = [];

for isite = sites_for_harmonics
    if exist(remotepath,'dir')
        [times{end+1},intensity] = grabdata(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
        times{end+1} = timestamp;
    end

    if log_trafo
        signals_raw{end+1} = log10(intensity);
    else
        signals_raw{end+1} = intensity;
    end
    
    signals{end+1} = signals_raw{end} - repmat(nanmean(signals_raw{end},2),1,size(signals_raw{end},2));
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

% DRUG DATA

dataPath = '2D_dose_response_drugsVSEGF_130903';

% site 20 (row 2 - col 1) is missing ...
sites_drug = [1:19 21:70];
lig_name = 'EGF';

times_drug = cell(0);
signals_drug = cell(0);
signals_drug_raw = cell(0);
celltype_drug = [];

for isite = sites_drug
    if exist(remotepath,'dir')
        [times_drug{end+1},intensity] = grabdata_drug(isite,dataPath);
    else
        load(['./Workspaces/site_' num2str(isite) '_only_' lig_name])
        times_drug{end+1} = timestamp;
    end

    if log_trafo
        signals_drug_raw{end+1} = log10(intensity);
    else
        signals_drug_raw{end+1} = intensity;
    end
    
    % Gate the signals at this place instead of at each call of grabdata_drug(isite,dataPath)
    signals_drug_raw{end} = signals_drug_raw{end}(:,(min(signals_drug_raw{end}) >= lower_gate_all) & (min(signals_drug_raw{end}(30:end,:)) >= lower_gate_late));
    
    signals_drug{end+1} = signals_drug_raw{end} - repmat(nanmean(signals_drug_raw{end},2),1,size(signals_drug_raw{end},2));
    
    celltype_drug = [celltype_drug ones(1,size(signals_drug{end},2))*isite];
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

nbasis = 40;

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

[tmp range_ind_min] = min(abs(timestamp_drug - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp_drug - time_range(2)));
range_ind_drug = range_ind_min:range_ind_max;

basis = create_bspline_basis([min([timestamp(range_ind(1)) timestamp_drug(range_ind_drug(1))]) max([timestamp(range_ind(end)) timestamp_drug(range_ind_drug(end))])], nbasis);
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

nharm = 3;
% c_signal_pcastr = pca_fd(smoothed_data, nharm);
c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!
% c_signal_pcastr = varmx_pca(c_signal_pcastr);

return

plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Plot: Eigenfunctions
close all

rowstocols = 0.5;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

flipharm = ones(1,nharm);
% flipharm(1:8) = [1 -1 1 -1 -1 1 1 1];

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

%% Remove all nharm from signal --> only stochastic oscillations remain
close all

evaluated_fd = eval_fd(c_signal_pcastr.fdhatfd,timestamp_drug(range_ind_drug));
c_signal_woNharm = c_signal_drug(range_ind_drug,:)-evaluated_fd(:,(size(c_signal,2)+1):end);

plot_sites = sites_drug;

rowstocols = 0.45;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

figure

for ip1 = 1:length(plot_sites)
    ip = plot_sites(ip1);
    subplot(nrows,ncols,ip)
    
    c_signal_single = c_signal_woNharm(:,celltype_drug == ip);
    
    first_n = 100;
    
    first_n = min(first_n,size(c_signal_single,2));
    
    plot(timestamp_drug(range_ind_drug),c_signal_single(:,1:first_n))
    prop_ind = (ip == plot_sites);
    title([lig_name ' ' num2str(site_lig_dose_drug(prop_ind)) site_inh_name_drug{prop_ind} num2str(site_inh_dose_drug(prop_ind))])
    
    set(gca,'XLim',time_range)
    
    hold on
    plot(timestamp_drug(range_ind_drug),nanmean(c_signal_single,2),'--k')
    
    plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
end

%% Generate spline fits to data-sets given in sites_for_harmonics (for remaining variation)
close all

nbasis = 40;

% ind_harm = ismember(celltype,sites_for_harmonics);
% ind_fit = ~ind_harm;

smoothed_data_woNharm = smooth_basis(timestamp_drug(range_ind_drug),c_signal_woNharm,basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal_drug(1,:),2)))
hold on

plot(smoothed_data_woNharm)
plot(timestamp_drug(range_ind_drug),c_signal_woNharm,'o')

%% Plot: Histogramm of distance to origin
close all

rad_dist_thres = 0.03;

figure
hold on

rowstocols = 0.45;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

posFig = get(gcf,'Position');
% posFig(3) = posFig(3)/2;
% posFig(4) = posFig(4)*2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

for ip1 = 1:length(plot_sites)
    ip = plot_sites(ip1);
    subplot(nrows,ncols,subplotpos(ip))
    
    baredges = linspace(0,max(radial_dist)+.01,21);
    bar(baredges,histc(radial_dist(celltype_drug == ip),baredges));
    
    prop_ind = (ip == plot_sites);
    title([lig_name ' ' num2str(site_lig_dose_drug(prop_ind)) site_inh_name_drug{prop_ind} num2str(site_inh_dose_drug(prop_ind))])
    
    set(gca,'XLim',[0 max(radial_dist)+.01])
    if ip == length(plot_sites)
        xlabel('radial distance')
    end
    if ip == 1
        ylabel('absolute frequency')
    end
    
    hold on
    
%     plot([rad_dist_thres rad_dist_thres],get(gca,'YLim'),'--')
end

%% Plot: Histogramm of distance to origin (new - overlayed)
close all

% groups = {[1 2 4 10], [4 10 17 57 64]};
% resort = {[2 3 1 4], [1 4 2 3 5]};

figure

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

sites_remain = 1:length(sites_drug);
uni_lig = unique(site_lig_dose_drug);
uni_drug = unique(site_inh_ind_drug(sites_remain));
drug_min = min(site_inh_dose_drug(sites_remain));
drug_max = max(site_inh_dose_drug(sites_remain));

rowstocols = 0.3;
nrows = length(uni_drug);
ncols = length(uni_lig);

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
%         ylabel('relative frequency')
%         ylabel('estimated density')
        ylabel('median(radial distance)')
    end
    if idrug == length(uni_drug) && ilig == 1
%         xlabel('radial distance')
        xlabel('log10(inhibitor dose)')
    end
    
    hold on
    
    mygroup = intersect(tmpind,tmpind2);
    
    legend_names = {};
    legendstyles = [];
    color = jet(length(mygroup));
    
    mx = [];
    mf = [];
    for ip = 1:length(mygroup)
        
        mydist = radial_dist(celltype_drug == mygroup(ip));
        if ~isempty(mydist)
            
            prop_ind = mygroup(ip);
            
%             baredges = linspace(0,max(radial_dist)+.01,26);
%             barheight = histc(mydist,baredges)./sum(celltype_drug == mygroup(ip));
%             legendstyles = [legendstyles plot(baredges,barheight,'Color',color(ip,:))];

%             [f,xi] = ksdensity(mydist);
%             mf = [mf sum(xi.*f)];
%             legendstyles = [legendstyles plot(xi,f,'Color',color(ip,:))];

%             legend_names{end+1} = num2str(site_inh_dose_drug(prop_ind));
            
            mf = [mf median(mydist)];
            mx = [mx site_inh_dose_drug(prop_ind)];
        end

    end
    
    mx = log10(mx);
    mx(mx==-Inf) = -1;
    plot(mx,mf,'k-')
    plot([-.8 -.8],get(gca,'YLim'),'b--')
    
%     legend(sps(end),legendstyles,legend_names)
    
%     set(gca,'XLim',[0 max(radial_dist)])
%     set(gca,'YLim',[0 .6])
%     set(gca,'YLim',[0 400])

    end
end
