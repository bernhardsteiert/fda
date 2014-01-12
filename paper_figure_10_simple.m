close all

% sites_all = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:70];
sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70]; % Without FGF
sites_colored = [11:17 64:70];

pcs = [2 3];

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f1,xfac,yfac,fontsize)

hold on

for isite = sites_all
    scores = fPCA(isite);
    plot(scores(2,:),scores(3,:),'.','Color',[.7 .7 .7])
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
    
    plot(scores(2,:),scores(3,:),'.','Color',mycolor)
    plotEllipsis(scores(2,:),scores(3,:),mycolor,.5)
    
end

xlim = [-.16 .24];
set(gca,'XLim',xlim)

axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
arrow([-.13 -.1],[.2 -.1],'Width',.5,'Length',7)
set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])
arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

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
    if sprop.lig_dose == 100
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
    plotEllipsis(scores(2,:),scores(3,:),colmap(isite == highdoses,:),.5)
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