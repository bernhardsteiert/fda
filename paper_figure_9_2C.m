%% Fit additional data-sets on harmonics generated from initial data-sets

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)


% sites = [17 41 42 44 57];
% input_names = {'IGF','HGF-MEKi','HGF-AKTi','HGF','EPR'};
% sites_for_harmonics = [17 41 42 44 57];

sites = [4 10 17 37 44 57 64];
input_names = {'EGF','No Lig','IGF','HRG','HGF','EPR','BTC'};
sites_for_harmonics = sites;
% sites_for_harmonics = [1 2 4 10 17 57 64];
% all ligands highest dose


times = cell(0);
signals = cell(0);
signals_raw = cell(0);
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
        signals_raw{end+1} = log10(intensity);
    else
        signals_raw{end+1} = intensity;
    end
    
    signals{end+1} = signals_raw{end} - repmat(nanmean(signals_raw{end},2),1,size(signals_raw{end},2));
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);


%% Plot raw data of defined ligands
close all
figure

c_signal_raw = cell2mat(signals_raw);

plot_ligs = 4;
plot_name = {'EGF 100 ng/ml - original'};
nrows = 2;
ncols = 3;

iplot = 1;
nplot = 10;
    
subplot(nrows,ncols,iplot)
plot(repmat(timestamp,1,sum(celltype == plot_ligs(iplot))),c_signal_raw(:,celltype == plot_ligs(iplot)),'g','color',[0.7 0.7 0.7])
colored_ind = find(celltype == plot_ligs(iplot));
hold on
plot(repmat(timestamp,1,nplot),c_signal_raw(:,colored_ind(1:nplot)))
plot(timestamp,nanmean(c_signal_raw(:,celltype == plot_ligs(iplot)),2),'color','k','LineWidth',2)
title(plot_name{iplot})

ylim = [-1 1]*.025+.005;
if ~log_trafo
    ylim = 10.^ylim;
end
set(gca,'XLim',[200 510],'YLim',ylim)
xlabel('time [min]')
ylabel('log_{10} FOXO3a Cyt/Nuc ratio');

subplot(nrows,ncols,iplot+1)
[tmp tmp_signal tmp_range] = radial_dist(plot_ligs);
plot(repmat(timestamp(tmp_range),1,sum(celltype == plot_ligs(iplot))),tmp_signal(:,celltype == plot_ligs(iplot)),'g','color',[0.7 0.7 0.7])
colored_ind = 1:size(tmp_signal,2);
hold on
plot(repmat(timestamp(tmp_range),1,nplot),tmp_signal(:,colored_ind(1:nplot)))
plot(timestamp(tmp_range),nanmean(tmp_signal(:,celltype == plot_ligs(iplot)),2),'color','k','LineWidth',2)
title('w/o mean + harm 1-3')

ylim = [-1 1]*.025+.005;
if ~log_trafo
    ylim = 10.^ylim;
end
set(gca,'XLim',[200 510],'YLim',ylim)


sites = [4 10 17 37 44 57 64];
input_names = {'EGF','No Lig','IGF','HRG','HGF','EPR','BTC'};

dists = [];
celltypeharm = [];

clear radial_dist

dists_mat = [];

for isite = sites
    radial_dists = radial_dist(isite);
    
    dists = [dists radial_dists];
    
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
    dists_mat = padconcatenation(dists_mat,radial_dists,1);

end

subplot(nrows,ncols,3)

resort = [2 3 6 4 1 5 7];
boxplot(dists_mat(resort,:)')
ylabel('Pulsatory strength')
set(gca,'XTick',1:7,'XTickLabel',input_names(resort))
set(gca,'YLim',[0 .04])

% %% Plot: Histogramm of distance to origin (new - overlayed)
% close all
% 
% groups = {sites_for_harmonics};
% resort = {1:7};
% 
% figure
% hold on
% 
% rowstocols = 0.6;
% nrows = ceil(length(groups)^rowstocols);
% ncols = ceil(length(groups) / nrows);
% 
% radial_dist = dists;
% 
% for ig = 1:length(groups)
%     g = subplot(nrows,ncols,ig);
%     hold on
%     
%     mygroup = groups{ig};
%     
%     legend_names = {};
%     legendstyles = nan(1,length(mygroup));
%     color = lines(length(mygroup));
%     color(4,:) = 0;
%     for ip = 1:length(mygroup)
% 
%         baredges = linspace(0,max(radial_dist)+.01,11);
%         barheight = histc(radial_dist(celltypeharm == mygroup(ip)),baredges)./sum(celltypeharm == mygroup(ip));
%         legendstyles(ip) = plot(baredges,barheight,'Color',color(resort{ig}(ip),:));
% 
% %         [f,xi] = ksdensity(radial_dist(celltypeharm == mygroup(ip)));
% %         legendstyles(ip) = plot(xi,f,'Color',color(resort{ig}(ip),:));
% 
%         legend_names{end+1} = (input_names{mygroup(ip) == sites});
% 
%     end
%     
%     legend(g,legendstyles,legend_names)
%     
%     set(gca,'XLim',[0 max(radial_dist)+.01])
% %     set(gca,'YLim',[0 200])
% 
%     xlabel('radial distance')
% %     ylabel('relative frequency')
%     ylabel('estimated density')
% 
% end


%% Plot traces with gray level corrsesponding to radial distance
plot_sites = [10 17 64];
plot_name = {'No Lig','IGF','BTC'};

nrows = 2;
ncols = 3;

posFig = get(gcf,'Position');
% posFig(3) = posFig(3)/2;
posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

ncolor = 201;
color = repmat(linspace(0,1,ncolor),3,1);

color = color(:,end:-1:1); % Gray scale - darkness depending on score

% TODO: As Pat suggested: Plot 4 representing lines / subplot in red - saturation depending on score

% color = (abs(1-color .* lines(ncolor)')).^2; % Colored - saturation depending on score
radial_space = linspace(min(dists),max(dists),ncolor);
[radial_dist_sorted ind_sort_radial] = sort(dists);

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip+3)
    hold on
    
    c_signal_single = c_signal_raw(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));
    radial_dist_single = radial_dist_sorted(celltypeharm(ind_sort_radial) == plot_sites(ip));
    
    for i = 1:size(c_signal_single,2)
        [tmp color_ind] = min(abs(radial_dist_single(i) - radial_space));
        plot(timestamp,c_signal_single(:,i),'Color',color(:,color_ind))
    end
    
    plot(get(gca,'XLim'),[0 0],'--k')
    
    title(plot_name{ip})
    
    set(gca,'XLim',[50 650])
    
    plot([120 120],[-0.04 0.04],'b--')
    
    if ip == 2
        xlabel('time [min]')
    end
    if ip == 1
        ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
    end
%     set(gca,'YLim',[-0.04 0.04])
end