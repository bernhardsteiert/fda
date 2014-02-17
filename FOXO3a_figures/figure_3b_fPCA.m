% Figure 3b: fPCA of early response in 184A1
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70]; % Without FGF

load('./Workspaces/scores_early')



highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);

    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
end

resort = [2 3 4 1 6 5];

highdoses = highdoses(resort);

figure;

plot(scores_early(2,~ismember(celltypes,highdoses)),scores_early(3,~ismember(celltypes,highdoses)),'o','MarkerSize',1,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none')
hold on

color_ind = 1;
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
legstr = cell(length(highdoses),1);
legh = [];
markers = {'o','s','v','d','^','>'};
for isite = highdoses([6 2 3 4 5 1])
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    
    scores = scores_early(:,celltypes == isite);
    legh(isite == highdoses) = plot(scores(2,:),scores(3,:),markers{isite == highdoses},'MarkerSize',5,'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerEdgeColor','none');
    plotEllipsis(scores(2,:),scores(3,:),colmap(isite == highdoses,:),.5);
end

xlim = [-.15 .25];
ylim = [-.05 .15];
set(gca,'XLim',xlim,'YLim',ylim)

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',1,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

ylabel('score of transient harmonic')
xlabel('score of sustained harmonic')

% set(gca,'CLim',[0 1])
% subplotpos = get(gca,'Position');
% colormap(colmap)
% colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0],'Position',[subplotpos(1)+subplotpos(3) subplotpos(2) .01 subplotpos(4)],'units','normalized') % Vertical colorbar
