% Figure 3b: fPCA of early response in 184A1
load('./Workspaces/scores_early_5basis_noFGF.mat')

sites_all = [4:10; 17:-1:11; 37:-1:31; 44:50; 57:-1:51; 64:69 10];
sites_high = [4 17 37 44 57 64];
sites_unst = [10 11 31 50 51];

myind = ~isnan(scores_early(1,:));

ps = 1:3;

for ipc = ps
    pvals_all = nan(size(sites_all));
    for i = 1:length(sites_all(:))
        isite = sites_all(i);
        pvals_all(i) = ranksum(scores_early(ipc,ismember(celltypes,sites_unst) & myind),scores_early(ipc,celltypes == isite & myind));
    end

    colmap = lines(size(sites_all,1));
    figure
    hold on
    legh = [];
    legstr = {};
    for isite = 1:size(sites_all,1)
        legh = [legh plot(pvals_all(isite,:),'Color',colmap(isite,:))];
        s = siteprop(sites_all(isite,1));
        legstr{end+1} = s.lig_name;
    end
    legend(legh,legstr,'Location','SouthWest')
    plot(get(gca,'XLim'),[1e-3 1e-3],'k--')
    set(gca,'YScale','log','XDir','reverse')
    title(sprintf('PC %i',ipc))
    xlabel('ligand dose')
    ylabel('p-value')
    set(gca,'XTick',1:7,'XTickLabel',[100 50 20 10 5 2.5 0])
end

return

med_lig = [];
iqr_lig = [];
pvals_lig = [];
for isite = 1:size(sites_all,1)
    pvals_lig = [pvals_lig ranksum(scores_early(ipc,myind),scores_early(ipc,ismember(celltypes,sites_all(isite,:)) & myind))];
    med_lig = [med_lig median(scores_early(ipc,ismember(celltypes,sites_all(isite,:)) & myind))];
    iqr_lig = [iqr_lig iqr(scores_early(ipc,ismember(celltypes,sites_all(isite,:)) & myind))/sqrt(sum(ismember(celltypes,sites_all(isite,:)) & myind))];
end
if ipc == 1
    % Compare against all
    med_lig = [med_lig median(scores_early(ipc,myind))];
    iqr_lig = [iqr_lig iqr(scores_early(ipc,myind))/sqrt(sum(myind))];
else
    % Compare against unstimulated
    med_lig = [med_lig median(scores_early(ipc,ismember(celltypes,sites_unst) & myind))];
    iqr_lig = [iqr_lig iqr(scores_early(ipc,ismember(celltypes,sites_unst) & myind))/sqrt(sum(ismember(celltypes,sites_unst) & myind))];
end
figure
errorbar(med_lig,iqr_lig)

pvals_high = [];
med_high = [];
iqr_high = [];
for isite = sites_high
    [tmp myp] = ttest2(scores_early(ipc,myind),scores_early(ipc,celltypes == isite & myind));
%     myp = ranksum(scores_early(ipc,myind),scores_early(ipc,celltypes == isite & myind));
    pvals_high = [pvals_high ];
    med_high = [med_high median(scores_early(ipc,celltypes == isite & myind))];
    iqr_high = [iqr_high iqr(scores_early(ipc,celltypes == isite & myind))/sqrt(sum(celltypes == isite & myind))];
end
med_high = [med_high med_lig(end)];
iqr_high = [iqr_high iqr_lig(end)];
figure
errorbar(med_high,iqr_high)
