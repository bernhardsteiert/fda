load('Workspaces/latePCsMedians.mat')
possible_doses = [0 2.5 5 10 20 50 100];

figure

for icol = 1:3
    subplot(1,3,icol)
    hold on
    legh = [];
    legstr = cell(length(highdoses),1);
    colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
    colmap = hsv2rgb(colmap(1:end-1,:));
    markers = {'o','s','v','d','^','>'};
    for irow = 1:length(highdoses)
        sprop = siteprop(highdoses(irow));
        legstr{irow} = sprop.lig_name(1:3);
        legh = [legh plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)];

        [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-nanmean(medians(1,:,icol)),1);
        plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + nanmean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)
        h=errorbar(0:length(possible_doses)-1,medians(:,irow,icol),medians(:,irow,icol+3),'Color',colmap(irow,:));
        set(h,'linestyle','none')
    end
    title(sprintf('Late PC%d',icol))
    set(gca,'XLim',[-.5 length(possible_doses)-.5])
    set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
    xlabel('Ligand dose [ng/ml]')
end