addpath('./Functions/')
load('./Workspaces/scores_early_5basis_noFGF_AKTi')
load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC_ATKi')
noInh = load('./Workspaces/scores_early_5basis_noFGF_newBTC');
noInhpuls = load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC');

pcs = [2 3];

sites_all = [17 37 44 4 64 57];
sites_akti = [19 39 42 2 62 59];

puls_thres = .3;

kswidth = 15;

mymin = min([scores_puls(:,1); noInhpuls.scores_puls(:,1)]);
mymax = max([scores_puls(:,1); noInhpuls.scores_puls(:,1)]);

highinds = ismember(celltypes,sites_akti);
[~,lambda,c] = boxcox([scores_puls(:,1); noInhpuls.scores_puls(:,1)]);
% --> c = 0 and lambda = 0.0500878494830929582581902081984
scores_puls_boxcox = boxcox_apply(scores_puls(highinds,1),lambda,c);
thres_trafo = boxcox_apply(puls_thres,lambda,c);
mymin_trafo = boxcox_apply(mymin,lambda,c);
mymax_trafo = boxcox_apply(mymax,lambda,c);
myrange = mymax_trafo - mymin_trafo;

scores_puls_boxcox = (scores_puls_boxcox-mymin_trafo)./myrange;
thres_trafo = (thres_trafo-mymin_trafo)./myrange;
% mymax_trafo = (mymax_trafo-mymin_trafo)./myrange;
% mymin_trafo = (mymin_trafo-mymin_trafo)./myrange;

figure

plot(noInh.scores_early(pcs(1),:),noInh.scores_early(pcs(2),:),'.','Color',[.7 .7 .7])
hold on

color_ind = 1;
colmap = hsv(length(sites_akti));
legstr = {};
for isite = sites_akti
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g %s',titstr,s.lig_dose,s.inh_name);
    legstr{end+1} = titstr;

    plot(scores_early(pcs(1),celltypes == isite),scores_early(pcs(2),celltypes == isite),'o','Color',colmap(isite == sites_akti,:),'MarkerFaceColor',colmap(isite == sites_akti,:))
    plotEllipsis(scores_early(pcs(1),celltypes == isite),scores_early(pcs(2),celltypes == isite),colmap(isite == sites_akti,:),.5);
end

% xlim = [-.4 .3];
% set(gca,'XLim',xlim)
% ylim = [-.08 .12];
% set(gca,'YLim',ylim)

% axisEqual(get(gcf,'Position'))

%     title(s.celltype)
ylabel(['PC ' num2str(pcs(2))])
% set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])

set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(sites_akti)),1-1./(2*length(sites_akti)),length(sites_akti)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar


figure
hold on

legh = [];
legstr = {};
for isite = sites_akti
    s = siteprop(isite);

    [f,xi] = ksdensity(scores_puls_boxcox(celltypes(highinds) == isite),'width',range(scores_puls(:,1))./kswidth);
    legh = [legh plot(xi,f,'Color',colmap(isite == sites_akti,:))];

    titstr = s.lig_name;
    titstr = sprintf('%s %g %s',titstr,s.lig_dose,s.inh_name);
    legstr{end+1} = titstr;

end
% set(gca,'XLim',[min(scores_puls_boxcox) max(scores_puls_boxcox)]+[-.2 .2]*range(scores_puls_boxcox))
%     set(gca,'XLim',xlim_all(icell,:))
%     title(s.celltype)
plot([thres_trafo thres_trafo],get(gca,'YLim'),'k--')
plot([0 0],get(gca,'YLim'),'k:')
plot([1 1],get(gca,'YLim'),'k:')
legend(legh,legstr)

xlabel('pulsatory score [au]')
ylabel('probability density')



kswidth = 40;

highinds = ismember(noInhpuls.celltypes,sites_all);
scores_puls_boxcox = boxcox_apply(noInhpuls.scores_puls(highinds,1),lambda,c);
scores_puls_boxcox = (scores_puls_boxcox-mymin_trafo)./myrange;

figure

plot(noInh.scores_early(pcs(1),:),noInh.scores_early(pcs(2),:),'.','Color',[.7 .7 .7])
hold on
plot(scores_early(pcs(1),:),scores_early(pcs(2),:),'.','Color',[.7 .7 .7])

color_ind = 1;
colmap = hsv(length(sites_all));
legstr = {};
for isite = sites_all
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    legstr{end+1} = titstr;

    plot(noInh.scores_early(pcs(1),noInh.celltypes == isite),noInh.scores_early(pcs(2),noInh.celltypes == isite),'o','Color',colmap(isite == sites_all,:),'MarkerFaceColor',colmap(isite == sites_all,:))
    plotEllipsis(noInh.scores_early(pcs(1),noInh.celltypes == isite),noInh.scores_early(pcs(2),noInh.celltypes == isite),colmap(isite == sites_all,:),.5);
end

% xlim = [-.4 .3];
% set(gca,'XLim',xlim)
% ylim = [-.08 .12];
% set(gca,'YLim',ylim)

% axisEqual(get(gcf,'Position'))

%     title(s.celltype)
ylabel(['PC ' num2str(pcs(2))])
% set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])

set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(sites_all)),1-1./(2*length(sites_all)),length(sites_all)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar


figure
hold on

legh = [];
legstr = {};
for isite = sites_all
    s = siteprop(isite);

    [f,xi] = ksdensity(scores_puls_boxcox(noInhpuls.celltypes(highinds) == isite),'width',range(noInhpuls.scores_puls(:,1))./kswidth);
    legh = [legh plot(xi,f,'Color',colmap(isite == sites_all,:))];

    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    legstr{end+1} = titstr;

end
% set(gca,'XLim',[min(scores_puls_boxcox) max(scores_puls_boxcox)]+[-.2 .2]*range(scores_puls_boxcox))
%     set(gca,'XLim',xlim_all(icell,:))
%     title(s.celltype)
plot([thres_trafo thres_trafo],get(gca,'YLim'),'k--')
plot([0 0],get(gca,'YLim'),'k:')
plot([1 1],get(gca,'YLim'),'k:')
legend(legh,legstr)

xlabel('pulsatory score [au]')
ylabel('probability density')
