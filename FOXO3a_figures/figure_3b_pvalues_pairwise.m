% Figure 3b: fPCA of early response in 184A1
load('./Workspaces/scores_early_5basis_noFGF_newBTC.mat')

sites_all = [17:-1:11; 37:-1:31; 44:50; 4:10; 64:69 10; 57:-1:51];
sites_high = [4 17 37 44 57 64];
sites_unst = [10 11 31 50 51];

myind = ~isnan(scores_early(1,:));

ps = 1:3;

colmap = hsv(size(sites_all,1));

for ipc = ps
    pvals_all = nan(length(sites_high));
    legstr = {};
    for i = 1:length(sites_high)
        isite = sites_high(i);
        for j = 1:length(sites_high)
            jsite = sites_high(j);
            pvals_all(i,j) = ranksum(scores_early(ipc,celltypes == jsite & myind),scores_early(ipc,celltypes == isite & myind));
        end
        s = siteprop(isite);
        legstr{end+1} = s.lig_name;
    end

    figure
    imagesc(log10(pvals_all))
    colorbar
    title(sprintf('PC %i',ipc))
    set(gca,'XTick',1:length(legstr),'XTickLabel',legstr)
    set(gca,'YTick',1:length(legstr),'YTickLabel',legstr)
    
    figure
    imagesc(log10(pvals_all)<-10)
    colormap('bone')
    title(sprintf('PC %i',ipc))
    set(gca,'XTick',1:length(legstr),'XTickLabel',legstr)
    set(gca,'YTick',1:length(legstr),'YTickLabel',legstr)
end
