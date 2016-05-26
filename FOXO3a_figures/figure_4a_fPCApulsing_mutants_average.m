addpath('./Functions/')
load('./Workspaces/scores_04-15_new')
close all;
extension = '04-15-2014';
pcs = [2 3];

highdoses_all = [1 24 25 48 49; 3 22 27 46 51; 5 20 29 44 53; 7 18 31 42 55; 9 16 33 40 57; 11 14 35 38 59];

kswidth_all = ones(1,6)*35;
% puls_thres = [.3 .5 .45 .2 .2 .2];
puls_thres = [.3 .5 .45 .1 .1 .1];

figure

plot(scores_all(pcs(1),:),scores_all(pcs(2),:),'.','Color',[.7 .7 .7])

markers = {'o','s','v','d','^','>'};
for icell = 1:size(highdoses_all,1)
    
    if icell == 4
        figure
        plot(scores_all(pcs(1),:),scores_all(pcs(2),:),'.','Color',[.7 .7 .7])
    end

    highdoses = highdoses_all(icell,:);

    hold on

    color_ind = 1;
    colmap = hsv(length(highdoses)+1);
    legstr = {};
    for isite = highdoses
        s = siteprop(isite,extension);
        titstr = s.lig_name;
        legstr{end+1} = titstr;
        
        if mod(icell,3) == 1
            reference(isite==highdoses,1:2) = [nanmean(scores_all(pcs(1),celltype == isite)),nanmean(scores_all(pcs(2),celltype == isite))];
            title(s.celltype(1:end-3))
            plot(nanmean(scores_all(pcs(1),celltype == isite)),nanmean(scores_all(pcs(2),celltype == isite)),markers{isite == highdoses},'Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerSize',12)
            plotEllipsis(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),colmap(isite == highdoses,:),2/sqrt(sum(~isnan(scores_all(pcs(1),celltype == isite)))));
        elseif mod(icell,3) == 0
            plot([reference(isite==highdoses,1) nanmean(scores_all(pcs(1),celltype == isite))],[reference(isite==highdoses,2) nanmean(scores_all(pcs(2),celltype == isite))],'--','Color','k')
            plot(nanmean(scores_all(pcs(1),celltype == isite)),nanmean(scores_all(pcs(2),celltype == isite)),markers{isite == highdoses},'Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerEdgeColor','w','MarkerSize',12)
            plotEllipsis(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),colmap(isite == highdoses,:),2/sqrt(sum(~isnan(scores_all(pcs(1),celltype == isite)))));
        else
            plot([reference(isite==highdoses,1) nanmean(scores_all(pcs(1),celltype == isite))],[reference(isite==highdoses,2) nanmean(scores_all(pcs(2),celltype == isite))],'-','Color',colmap(isite == highdoses,:),'LineWidth',2)
            plot(nanmean(scores_all(pcs(1),celltype == isite)),nanmean(scores_all(pcs(2),celltype == isite)),markers{isite == highdoses},'Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerEdgeColor','k','MarkerSize',12)
            plotEllipsis(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),colmap(isite == highdoses,:),2/sqrt(sum(~isnan(scores_all(pcs(1),celltype == isite)))));
        end
        
    end

    xlim = [-0.1 .2];
    set(gca,'XLim',xlim)
    ylim = [-.02 .07];
    set(gca,'YLim',ylim)

    ylabel(['fPC ' num2str(pcs(2))])
    xlabel(['fPC ' num2str(pcs(1))])

    set(gca,'CLim',[0 1])
    colormap(colmap(1:end-1,:))
    colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', 0) % Vertical colorbar
    
end

figure


count = 1;
errorb = nan(length(highdoses_all(:)),4);
mutstr = {};
for icell = 1:size(highdoses_all,1)
    
    rng(0) % Make sure that bootstrap samples are reproducible
    ratio_fun = @(x) sum(x > puls_thres(icell)) / length(x);
    
    if icell == 4
        errorbar(errorb(:,1),errorb(:,2),errorb(:,2)-errorb(:,3),errorb(:,4)-errorb(:,2),'LineStyle','none','Color','k');
        figure
        errorb = nan(length(highdoses_all(:)),4);
    end
    hold on

    highdoses = highdoses_all(icell,:);
    highinds = ismember(celltype,highdoses);
    dists2 = dists(highinds);

    legh = [];
    legstr = {};
    
    for i = 1:length(highdoses)
        isite = highdoses(i);
        s = siteprop(isite,extension);
        
        errorb(count,1) = mod(icell-1,3)*size(highdoses_all,2)+i;
        errorb(count,2) = sum(dists2(celltype(highinds) == isite) > puls_thres(icell))./sum(celltype(highinds) == isite);
        
        legh = [legh bar(errorb(count,1),errorb(count,2),'FaceColor',colmap(isite == highdoses,:))];
        legstr{end+1} = s.lig_name;
        
        errorb(count,3:4) = bootci(2000,{ratio_fun,dists2(celltype(highinds) == isite)},'alpha',.32);
        
        count = count + 1;
    end
    mutstr{end+1} = s.celltype(8:end);
    set(gca,'XTick',3:5:13,'XTickLabel',mutstr);

    title(s.celltype(1:end-8))
    legend(legh,legstr)

%     xlabel('ligand')
    ylabel('fraction of pulsing cells')

end

errorbar(errorb(:,1),errorb(:,2),errorb(:,2)-errorb(:,3),errorb(:,4)-errorb(:,2),'LineStyle','none','Color','k');