%% Make analysis of pulsatory behavior different ligand concentrations

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

log_trafo = 1; % log-transform signal
time_range = [200 650];
time_range = [200 510];

% Get properties of sites by calling siteprop(site)
sites_for_harmonics = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:70];
% Note: Site 70 is blank and therefore a copy of site 30

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


timestamp = times{1}; % same time sampling for all data sets
c_signal_raw = cell2mat(signals_raw);
c_signal = cell2mat(signals);

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


%% Generate spline fits to data-sets given in sites_for_harmonics
% close all

nbasis = 40;

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),basis);

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

evaluated_fd = eval_fd(c_signal_pcastr.fdhatfd,timestamp(range_ind));
c_signal_woNharm = c_signal(range_ind,:)-evaluated_fd;

plot_sites = sites_for_harmonics;

rowstocols = 0.5;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

figure

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip)
    
    c_signal_single = c_signal_woNharm(:,celltype == plot_sites(ip));
    
    first_n = 100;
    
    first_n = min(first_n,size(c_signal_single,2));
    
    plot(timestamp(range_ind),c_signal_single(:,1:first_n))
    title([site_lig_name{ip} num2str(site_lig_dose(ip))])
    
    set(gca,'XLim',time_range)
    
    hold on
    plot(timestamp(range_ind),nanmean(c_signal_single,2),'--k')
    
    plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
end

%% Generate spline fits to data-sets given in sites_for_harmonics (for remaining variation)
close all

nbasis = 40;

% ind_harm = ismember(celltype,sites_for_harmonics);
% ind_fit = ~ind_harm;

smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_woNharm,basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,:),2)))
hold on

plot(smoothed_data_woNharm)
plot(timestamp(range_ind),c_signal_woNharm,'o')

%% Plot: Histogramm of distance to origin
close all

rad_dist_thres = 0.03;

figure
hold on

rowstocols = 0.5;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

posFig = get(gcf,'Position');
% posFig(3) = posFig(3)/2;
% posFig(4) = posFig(4)*2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip)
    
    baredges = linspace(0,max(radial_dist)+.01,21);
    bar(baredges,histc(radial_dist(celltype == plot_sites(ip)),baredges));
    
    title([site_lig_name{ip} num2str(site_lig_dose(ip))])
    
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

figure

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

sites_remain = sites_for_harmonics;
uni_lig = unique(site_lig_ind);

rowstocols = 0.5;
nrows = ceil(length(uni_lig)^rowstocols);
ncols = ceil(length(uni_lig) / nrows);

posFig = get(gcf,'Position');
% posFig(4) = posFig(4)/3;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./20);
sps = [];

for ilig = uni_lig

    tmp = subplot(nrows,ncols,ilig);
    sps = [sps tmp];
    hold on
    tmpind = find(site_lig_ind(1:length(sites_remain)) == ilig);
    title(site_lig_name{tmpind(1)})

    if ilig == 1
%         ylabel('relative frequency')
%         ylabel('estimated density')
        ylabel('median(pulsatory strength)')
    end
    if ilig == uni_lig(end)
%         xlabel('radial distance')
        xlabel('log10(ligand dose)')
    end

    hold on

    legend_names = {};
    legendstyles = [];
    color = jet(length(tmpind));

    mx = [];
    mf = [];
    for ip = 1:length(tmpind)

        mydist = radial_dist(celltype == sites_for_harmonics(tmpind(ip)));
        if ~isempty(mydist)

            prop_ind = tmpind(ip);

            baredges = linspace(0,max(radial_dist)+.01,26);
            barheight = histc(mydist,baredges)./sum(celltype == sites_for_harmonics(tmpind(ip)));
%             legendstyles = [legendstyles plot(baredges,barheight,'Color',color(ip,:))];

%             [f,xi] = ksdensity(mydist);
%             legendstyles = [legendstyles plot(xi,f,'Color',color(ip,:))];

%             legend_names{end+1} = num2str(site_lig_dose(prop_ind));

            mf = [mf median(mydist)];
            mx = [mx site_lig_dose(prop_ind)];
        end

    end

    mx = log10(mx);
    mx(mx==-Inf) = 0;
    plot(mx,mf,'x')
    set(gca,'XLim',[-.1 2.1],'YLim',[.005 .02])
    plot([0.2 .2],get(gca,'YLim'),'b--')
    
    [axb s] = polyfit(mx(2:end),mf(2:end),1);
    plot(mx,mx*axb(1) + axb(2),'k--')
    
    Rinv = inv(s.R);
    covmat = (Rinv*Rinv')*s.normr^2/s.df;
    plot(mx,mx*(axb(1)+sqrt(covmat(1))) + axb(2),'k:')
    plot(mx,mx*(axb(1)-sqrt(covmat(1))) + axb(2),'k:')

%     legend(sps(end),legendstyles,legend_names)

%     set(gca,'XLim',[0 max(radial_dist)])
%     set(gca,'YLim',[0 .6])
%     set(gca,'YLim',[0 120])

end


%% Plot: Scatterplot to see if radial_dist and inital values are correlated
close all

initial_timerange = [50 115];
[tmp initial_range_ind_min] = min(abs(timestamp - initial_timerange(1)));
[tmp initial_range_ind_max] = min(abs(timestamp - initial_timerange(2)));
initial_range_ind = initial_range_ind_min:initial_range_ind_max;

figure

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

sites_remain = sites_for_harmonics;
uni_lig = unique(site_lig_ind);

rowstocols = 0.5;
nrows = ceil(length(uni_lig)^rowstocols);
ncols = ceil(length(uni_lig) / nrows);

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/1.5;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);
sps = [];

for ilig = uni_lig

    tmp = subplot(nrows,ncols,ilig);
    sps = [sps tmp];
    hold on
    tmpind = find(site_lig_ind(1:length(sites_remain)) == ilig);
    title(site_lig_name{tmpind(1)})

    if ilig == 1
        ylabel('initial value')
    end
    if ilig == uni_lig(end)
        xlabel('radial distance')
    end

    hold on

    legend_names = {};
    legendstyles = [];
    color = jet(length(tmpind));

    mx = [];
    mf = [];
    for ip = 1:length(tmpind)

        mydist = radial_dist(celltype == sites_for_harmonics(tmpind(ip)));
        if ~isempty(mydist)

            prop_ind = tmpind(ip);
            legendstyles = [legendstyles plot(mydist,mean(c_signal_raw(initial_range_ind,celltype == sites_for_harmonics(tmpind(ip))),1),'.','Color',color(ip,:))];

            legend_names{end+1} = num2str(site_lig_dose(prop_ind));
            
        end

    end

    legend(sps(end),legendstyles,legend_names)
    
    mean_init = mean(c_signal_raw(initial_range_ind,:));
    set(gca,'XLim',[0 max(radial_dist)])
    set(gca,'YLim',[min(mean_init) max(mean_init)])

end

%% Plot: Boxplot
close all

alpha = 0.001;

figure

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

sites_remain = sites_for_harmonics;
uni_lig = unique(site_lig_ind);

rowstocols = 0.5;
nrows = ceil(length(uni_lig)^rowstocols);
ncols = ceil(length(uni_lig) / nrows);

posFig = get(gcf,'Position');
% posFig(4) = posFig(4)/3;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./20);
sps = [];

for ilig = uni_lig

    tmp = subplot(nrows,ncols,ilig);
    sps = [sps tmp];
    hold on
    tmpind = find(site_lig_ind(1:length(sites_remain)) == ilig);
    title(site_lig_name{tmpind(1)})

%     if ilig == 1
        ylabel('radial distance')
%     end
    if ilig == uni_lig(end)
        xlabel('ligand dose')
    end

    hold on

    legend_names = {};
    legendstyles = [];
    color = jet(length(tmpind));

    mf = [];
    for ip = length(tmpind):-1:1
        
        mydist = radial_dist(celltype == sites_for_harmonics(tmpind(ip)));
        if ~isempty(mydist)
            if length(mydist) > size(mf,1)
                mf = [mf; nan(length(mydist)-size(mf,1),size(mf,2))];
            else
                mydist = [mydist nan(1,size(mf,1)-length(mydist))];
            end
            mf = [mf mydist'];
            dist1 = mf(:,1);
            % Wilcoxon Test:
            if ranksum(dist1(~isnan(dist1)),mydist(~isnan(mydist))') < alpha
                plot(size(mf,2),.055,'k*','LineWidth',1)
            end
            
        end
        
        legend_names{end+1} = num2str(site_lig_dose(tmpind(ip)));

    end
    
    boxplot(mf)
    set(gca,'XTick',1:7,'XTickLabel',legend_names)

    set(gca,'YLim',[0 max(radial_dist)])

end

subplot(nrows,ncols,ilig+1);
set(gca,'Visible','off')
text(0,1,'* p < 0.001 against NS')