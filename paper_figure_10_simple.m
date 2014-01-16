close all

% sites_all = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:70];
sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70]; % Without FGF
sites_colored = [11:17 64:70];

pcs = [1 2];

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f1,xfac,yfac,fontsize)

hold on

for isite = sites_all
    scores = fPCA(isite);
    plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
end

ncolor = 201;
colmap = jet(ncolor);
dose_range = [0 100];
dose_inds = 10.^linspace(log10(max([dose_range(1) 1])),log10(dose_range(2)),floor(ncolor/2));

firstprop = siteprop(sites_colored(1));
for isite = sites_colored
    scores = fPCA(isite);
    sprop = siteprop(isite);
    
    colfac = 2*(sprop.lig_index == firstprop.lig_index)-1;
    [tmp colind] = min(abs(sprop.lig_dose - dose_inds));
    mycolor = colmap(ceil(ncolor/2) + colfac*colind,:);
    
    plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',mycolor)
    plotEllipsis(scores(pcs(1),:),scores(pcs(2),:),mycolor,.5);
    
end

% xlim = [-.16 .24];
% set(gca,'XLim',xlim)

% axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
% arrow([-.13 -.1],[.2 -.1],'Width',.5,'Length',7)
set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])
% arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

clim = log10(dose_range);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(colmap)
colorbar('YTick',log10([1 10 100]),'YTickLabel',{'BTC','No Stim','IGF'}) % Vertical colorbar

% -------------------------------------------------------------------------

f2 = figure;
hold on

xfac = 1;
yfac = .6;

setFigure(f2,xfac,yfac,fontsize)

highdoses = [];
for isite = sites_all
    scores = fPCA(isite);
    plot(scores(2,:),scores(3,:),'.','Color',[.7 .7 .7])
    
    sprop = siteprop(isite);
    if sprop.lig_dose == 50
        highdoses = [highdoses isite];
    end
end

resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

color_ind = 1;
colmap = hsv(length(highdoses));
legstr = {};
for isite = highdoses
    sprop = siteprop(isite);
    legstr{end+1} = sprop.lig_name(1:3);
    
    scores = fPCA(isite);
    plot(scores(2,:),scores(3,:),'.','Color',colmap(isite == highdoses,:))
    plotEllipsis(scores(2,:),scores(3,:),colmap(isite == highdoses,:),.5);
end

xlim = [-.16 .24];
set(gca,'XLim',xlim)

axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
% arrow([-.13 -.1],[.2 -.1],'Width',.5,'Length',7)
set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])
% arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar

% -------------------------------------------------------------------------

f3 = figure;

xfac = 1.5;
yfac = 1;
fontsize = 10;

setFigure(f3,xfac,yfac,fontsize)

doses = [];
heterogeneity_early = [];
heterogeneity_late = [];
leg_str = {};

for isite = sites_all
    sprop = siteprop(isite);
    if size(heterogeneity_late,1) < sprop.lig_index
        heterogeneity_late = [heterogeneity_late; nan(size(heterogeneity_late,1)-sprop.lig_index+1,size(heterogeneity_late,2))];
        leg_str{sprop.lig_index} = sprop.lig_name;
    end
    
    ind_dose = find(ismember(doses,sprop.lig_dose));
    if isempty(ind_dose)
        heterogeneity_late = [heterogeneity_late nan(size(heterogeneity_late,1),1)];
        ind_dose = size(heterogeneity_late,2);
        doses = [doses sprop.lig_dose];
    end
    
    scores = fPCA(isite);
    scores_medsub = scores - repmat(median(scores,2),1,size(scores,2));
    [tmp indmedsort] = sort(scores_medsub(2,:).^2 + scores_medsub(3,:).^2);
    scores = scores(:,indmedsort(1:round(length(indmedsort)*50/68))); % kind-of-robust
    
    heterogeneity_early(sprop.lig_index,ind_dose) = ellipsis_area(scores(2,:),scores(3,:));
    heterogeneity_late(sprop.lig_index,ind_dose) = iqr(radial_dist(isite));

end

lig_inds = 1:size(heterogeneity_late,1);
lig_inds = lig_inds(sum(heterogeneity_late,2) > 0);
heterogeneity_early = heterogeneity_early(sum(heterogeneity_late,2) > 0,:);
heterogeneity_late = heterogeneity_late(sum(heterogeneity_late,2) > 0,:);

resort = [2 3 4 1 6 5];
heterogeneity_early = heterogeneity_early(resort,end:-1:1);
heterogeneity_late = heterogeneity_late(resort,end:-1:1);
doses = doses(end:-1:1);
lig_inds = lig_inds(resort);

subplot(1,2,1)
imagesc(heterogeneity_early)

set(gca,'XTick',1:size(heterogeneity_late,2),'XTickLabel',doses)
set(gca,'YTick',1:size(heterogeneity_late,1),'YTickLabel',leg_str(lig_inds))

colormap('hot')
colorbar

xlabel('Ligand dose')

title('Heterogeneity (early)')

subplot(1,2,2)
imagesc(heterogeneity_late)

set(gca,'XTick',1:size(heterogeneity_late,2),'XTickLabel',doses)
set(gca,'YTick',1:size(heterogeneity_late,1),'YTickLabel',leg_str(lig_inds))

colormap('hot')
colorbar

xlabel('Ligand dose')

title('Heterogeneity (late)')
        