load('Workspaces/scores_03242014')
addpath('./Functions/')

sites_all = [1:60];

nrows = 6;
ncols = 10;

% Sort by ratio of pulsing cells
puls_thres = .8;

labels = {};
ratio_puls = nan(size(sites_all));
ratio_puls_mat = nan(ncols,nrows);

mean_amp = nan(size(sites_all));
mean_amp_mat = nan(ncols,nrows);

mean_peakdur = nan(size(sites_all));
mean_peakdur_mat = nan(ncols,nrows);

for isite = sites_all
    ratio_puls(isite) = sum(dists(celltype == isite) > puls_thres) / sum(~isnan(dists(celltype == isite)));
    ratio_puls_mat(subplotpos(isite,ncols)) = ratio_puls(isite);
end

% Pulsing(Ligand dose)
close all
imagesc(ratio_puls_mat')
colorbar
possible_doses = [100 100/2.5 100/2.5^2 100/2.5^3 100/2.5^4 100/2.5^5 100/2.5^6 100/2.5^7 100/2.5^8 0];
possible_ligands = {'IGF','HRG','HGF','EGF','BTC','EPR'};
set(gca,'XTick',1:length(possible_doses),'XTickLabel',possible_doses)
set(gca,'YTick',1:length(possible_ligands),'YTickLabel',possible_ligands)
colormap(cbrewer('seq','YlGnBu',201))

possible_doses_main = [100 50 20 10 5 2.5 0];
log_doses = log10(possible_doses);
log_doses(end) = log10(100/2.5^9);
log_doses_main = log10(possible_doses_main);
log_doses_main(end) = log_doses(end);

ratio_puls_mat_interp = nan(size(ratio_puls_mat,2),length(log_doses_main));
for ilig = 1:size(ratio_puls_mat,2)
    ratio_puls_mat_interp(ilig,:) = interp1(log_doses,ratio_puls_mat(:,ilig)',log_doses_main);
end

figure
imagesc(ratio_puls_mat_interp)
colorbar
set(gca,'XTick',1:length(possible_doses_main),'XTickLabel',possible_doses_main)
set(gca,'YTick',1:length(possible_ligands),'YTickLabel',possible_ligands)

colormap(cbrewer('seq','YlGnBu',201))