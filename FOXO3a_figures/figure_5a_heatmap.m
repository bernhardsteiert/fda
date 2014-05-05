% Figure 5a: Heatmap EGF vs. MEKi
addpath('./Functions/')
load('Workspaces/dists_04182014')
extension = '04-18-2014';
sites_all = [1:39 41:72];

nrows = 6;
ncols = 12;

figure

puls_thres = .55;

ratio_puls = nan(size(sites_all));
ratio_puls_mat = nan(ncols,nrows);

mean_amp = nan(size(sites_all));
mean_amp_mat = nan(ncols,nrows);

mean_peakdur = nan(size(sites_all));
mean_peakdur_mat = nan(ncols,nrows);

for isite = sites_all
    ratio_puls(isite) = sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite);
    
    % Inside pulsing subpopulation
    mean_amp(isite) = mean(amps(dists > puls_thres & celltype == isite));
    mean_peakdur(isite) = mean(peakdurs(dists > puls_thres & celltype == isite));
    
%     % Complete population
%     mean_amp(isite) = nanmean(amps(celltype == isite));
%     mean_peakdur(isite) = nanmean(peakdurs(celltype == isite));
    ratio_puls_mat(subplotpos(isite,ncols)) = ratio_puls(isite);
    mean_amp_mat(subplotpos(isite,ncols)) = mean_amp(isite);
    mean_peakdur_mat(subplotpos(isite,ncols)) = mean_peakdur(isite);
end

[ratio_puls_sorted ind_puls_sorted] = sort(ratio_puls,'descend');

[X,Y] = meshgrid(1:size(ratio_puls_mat,2), 1:1:size(ratio_puls_mat,1)); 

valid = ~isnan(ratio_puls_mat); 
M = griddata(X(valid),Y(valid),ratio_puls_mat(valid),X,Y);
zlab = 'Ratio of pulsing cells';

colmap = cbrewer('seq','YlGnBu',201);

figure

imagesc(M(1:6,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['MCF10A - ' zlab])
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
colorbar

figure

imagesc(M(12:-1:7,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['184A1 - ' zlab])
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
colorbar

valid = ~isnan(mean_amp_mat); 
M = griddata(X(valid),Y(valid),mean_amp_mat(valid),X,Y); 
zlab = 'Amplitude';

figure

imagesc(M(1:6,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['MCF10A - ' zlab])
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
colorbar

figure

imagesc(M(12:-1:7,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['184A1 - ' zlab])
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
colorbar

valid = ~isnan(mean_peakdur_mat); 
M = griddata(X(valid),Y(valid),mean_peakdur_mat(valid),X,Y); 
zlab = 'Peak duration';

imagesc(M(1:6,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['MCF10A - ' zlab])
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
colorbar

figure

imagesc(M(12:-1:7,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['184A1 - ' zlab])
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
colorbar