addpath('./Functions/')
load('./Workspaces/scores_early_5basis_noFGF_AKTi')
meki = load('./Workspaces/scores_early_5basis_noFGF_MEKi');
mekipuls = load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC_MEKi');
load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC_ATKi')
noInh = load('./Workspaces/scores_early_5basis_noFGF_newBTC');
noInhpuls = load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC');
pcs = [2 3];

sites_all = [17 37 44 4 64];
sites_akti = [19 39 42 2 62];
sites_meki = [20 40 41 1 61];

puls_thres = .3;

highinds = ismember(celltypes,sites_akti);
scores_puls = scores_puls(highinds,1);
highindsNoInh = ismember(noInhpuls.celltypes,sites_all);
noInhpuls.scores_puls = noInhpuls.scores_puls(highindsNoInh,1);
highindsMEKi = ismember(mekipuls.celltypes,sites_meki);
mekipuls.scores_puls = mekipuls.scores_puls(highindsMEKi,1);

figure

% plot(noInh.scores_early(pcs(1),:),noInh.scores_early(pcs(2),:),'.','Color',[.7 .7 .7])
% plot(noInh.scores_early(pcs(1),highindsNoInh),noInh.scores_early(pcs(2),highindsNoInh),'.','Color',[.7 .7 .7]) % WT
hold on
% plot(scores_early(pcs(1),:),scores_early(pcs(2),:),'.','Color',[.7 .7 .7])
% plot(scores_early(pcs(1),highinds),scores_early(pcs(2),highinds),'.','Color',[.7 .7 .7]) % AKTi
% plot(meki.scores_early(pcs(1),highindsMEKi),meki.scores_early(pcs(2),highindsMEKi),'.','Color',[.7 .7 .7]) % MEKi

color_ind = 1;
colmap = hsv(length(sites_akti)+1);
legstr = {};
markers = {'o','s','v','d','^','>'};
for i = 1:length(sites_akti)
    isite = sites_akti(i);
    s = siteprop(isite);
    titstr = s.lig_name;
%     titstr = sprintf('%s %g %s',titstr,s.lig_dose,s.inh_name);
    legstr{end+1} = titstr;
    
    plot(nanmean(scores_early(pcs(1),celltypes == isite)),nanmean(scores_early(pcs(2),celltypes == isite)),markers{isite == sites_akti},'Color',colmap(isite == sites_akti,:),'MarkerFaceColor',colmap(isite == sites_akti,:),'MarkerEdgeColor','w','MarkerSize',12)
    plotEllipsis(scores_early(pcs(1),celltypes == isite),scores_early(pcs(2),celltypes == isite),colmap(isite == sites_akti,:),2/sqrt(sum(~isnan(scores_early(pcs(1),celltypes == isite)))));
    
    isite3 = sites_meki(i);
    plot(nanmean(meki.scores_early(pcs(1),meki.celltypes == isite3)),nanmean(meki.scores_early(pcs(2),meki.celltypes == isite3)),markers{isite3 == sites_meki},'Color',colmap(isite3 == sites_meki,:),'MarkerFaceColor',colmap(isite3 == sites_meki,:),'MarkerEdgeColor','k','MarkerSize',12)
    plotEllipsis(meki.scores_early(pcs(1),meki.celltypes == isite3),meki.scores_early(pcs(2),meki.celltypes == isite3),colmap(isite3 == sites_meki,:),2/sqrt(sum(~isnan(meki.scores_early(pcs(1),meki.celltypes == isite3)))));
    isite2 = sites_all(i);
    plot(nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite2)),nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite2)),markers{isite2 == sites_all},'Color',colmap(isite2 == sites_all,:),'MarkerFaceColor',colmap(isite2 == sites_all,:),'MarkerSize',12)
    plotEllipsis(noInh.scores_early(pcs(1),noInh.celltypes == isite2),noInh.scores_early(pcs(2),noInh.celltypes == isite2),colmap(isite2 == sites_all,:),2/sqrt(sum(~isnan(noInh.scores_early(pcs(1),noInh.celltypes == isite2)))));
    
    plot([nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite2)) nanmean(meki.scores_early(pcs(1),meki.celltypes == isite3))],[nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite2)) nanmean(meki.scores_early(pcs(2),meki.celltypes == isite3))],'-','Color',colmap(isite == sites_akti,:),'LineWidth',2)
    plot([nanmean(scores_early(pcs(1),celltypes == isite)) nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite2))],[nanmean(scores_early(pcs(2),celltypes == isite)) nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite2))],'k--','Color','k')
%     plotEllipsis(scores_early(pcs(1),celltypes == isite),scores_early(pcs(2),celltypes == isite),colmap(isite == sites_akti,:),.5);
end
set(gca,'XLim',[-0.3 0.25],'YLim',[-0.02 0.1])
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
colormap(colmap(1:end-1,:))
colorbar('YTick',linspace(1./(2*length(sites_akti)),1-1./(2*length(sites_akti)),length(sites_akti)),'YTickLabel',legstr,'TickLength', 0) % Vertical colorbar
% return

figure
hold on

rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x) sum(x > puls_thres) / length(x);
errorb = nan(3*length(sites_akti),4);

legh = [];
legstr = {};
for i = 1:length(sites_akti)
    isite = sites_akti(i);
    s = siteprop(isite);
    
    iinh = 3;
    errorb((i-1)*3+iinh,1) = mod(iinh-1,3)*(length(sites_akti)+1)+i;
    errorb((i-1)*3+iinh,2) = sum(scores_puls(celltypes(highinds) == isite) > puls_thres)./sum(celltypes(highinds) == isite);
    errorb((i-1)*3+iinh,3:4) = bootci(2000,{ratio_fun,scores_puls(celltypes(highinds) == isite)},'alpha',.32);
    legh = [legh bar(errorb((i-1)*3+iinh,1),errorb((i-1)*3+iinh,2),'FaceColor',colmap(isite == sites_akti,:))];
    iinh = 1;
    isiteH = sites_all(i);
    errorb((i-1)*3+iinh,1) = mod(iinh-1,3)*(length(sites_akti)+1)+i;
    errorb((i-1)*3+iinh,2) = sum(noInhpuls.scores_puls(noInhpuls.celltypes(highindsNoInh) == isiteH) > puls_thres)./sum(noInhpuls.celltypes(highindsNoInh) == isiteH);
    errorb((i-1)*3+iinh,3:4) = bootci(2000,{ratio_fun,noInhpuls.scores_puls(noInhpuls.celltypes(highindsNoInh) == isiteH)},'alpha',.32);
    bar(errorb((i-1)*3+iinh,1),errorb((i-1)*3+iinh,2),'FaceColor',colmap(isiteH == sites_all,:));
    iinh = 2;
    isiteM = sites_meki(i);
    errorb((i-1)*3+iinh,1) = mod(iinh-1,3)*(length(sites_akti)+1)+i;
    errorb((i-1)*3+iinh,2) = sum(mekipuls.scores_puls(mekipuls.celltypes(highindsMEKi) == isiteM) > puls_thres)./sum(mekipuls.celltypes(highindsMEKi) == isiteM);
    errorb((i-1)*3+iinh,3:4) = bootci(2000,{ratio_fun,mekipuls.scores_puls(mekipuls.celltypes(highindsMEKi) == isiteM)},'alpha',.32);
    bar(errorb((i-1)*3+iinh,1),errorb((i-1)*3+iinh,2),'FaceColor',colmap(isiteM == sites_meki,:));
    legstr{end+1} = s.lig_name;
%     [f,xi] = ksdensity(scores_puls(celltypes(highinds) == isite));
%     legh = [legh plot(xi,f,'Color',colmap(isite == sites_akti,:))];

%     titstr = s.lig_name;
%     titstr = sprintf('%s %g %s',titstr,s.lig_dose,s.inh_name);

end
%
%     set(gca,'XLim',xlim_all(icell,:))
title('184A1')
set(gca,'XTick',3:5:13,'XTickLabel',{'WT','MEKi','AKTi'});
% plot([thres_trafo thres_trafo],get(gca,'YLim'),'k--')
% plot([0 0],get(gca,'YLim'),'k:')
% plot([1 1],get(gca,'YLim'),'k:')
legend(legh,legstr)

% xlabel('ligand')
ylabel('fraction of pulsing cells')

errorbar(errorb(:,1),errorb(:,2),errorb(:,2)-errorb(:,3),errorb(:,4)-errorb(:,2),'LineStyle','none','Color','k');

return
% 
% 
% kswidth = 40;
% 
% highinds = ismember(noInhpuls.celltypes,sites_all);
% scores_puls_boxcox = boxcox_apply(noInhpuls.scores_puls(highinds,1),lambda,c);
% scores_puls_boxcox = (scores_puls_boxcox-mymin_trafo)./myrange;

% figure

% plot(noInh.scores_early(pcs(1),:),noInh.scores_early(pcs(2),:),'.','Color',[.7 .7 .7])
% hold on
% plot(scores_early(pcs(1),:),scores_early(pcs(2),:),'.','Color',[.7 .7 .7])

color_ind = 1;
colmap = hsv(length(sites_all));
legstr = {};
for isite = sites_all
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    legstr{end+1} = titstr;

    plot(nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite)),nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite)),markers{isite == sites_all},'Color',colmap(isite == sites_all,:),'MarkerFaceColor',colmap(isite == sites_all,:))
%     plotEllipsis(noInh.scores_early(pcs(1),noInh.celltypes == isite),noInh.scores_early(pcs(2),noInh.celltypes == isite),colmap(isite == sites_all,:),.5);
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


% figure
% hold on
% 
% legh = [];
% legstr = {};
% for isite = sites_all
%     s = siteprop(isite);
% 
%     [f,xi] = ksdensity(scores_puls_boxcox(noInhpuls.celltypes(highinds) == isite),'width',range(noInhpuls.scores_puls(:,1))./kswidth);
%     legh = [legh plot(xi,f,'Color',colmap(isite == sites_all,:))];
% 
%     titstr = s.lig_name;
%     titstr = sprintf('%s %g',titstr,s.lig_dose);
%     legstr{end+1} = titstr;
% 
% end
% % set(gca,'XLim',[min(scores_puls_boxcox) max(scores_puls_boxcox)]+[-.2 .2]*range(scores_puls_boxcox))
% %     set(gca,'XLim',xlim_all(icell,:))
% %     title(s.celltype)
% plot([thres_trafo thres_trafo],get(gca,'YLim'),'k--')
% plot([0 0],get(gca,'YLim'),'k:')
% plot([1 1],get(gca,'YLim'),'k:')
% legend(legh,legstr)
% 
% xlabel('pulsatory score [au]')
% ylabel('probability density')
