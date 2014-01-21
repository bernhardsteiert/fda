close all

load('Workspaces/site_2_12-08-2013.mat')
timeshift = -100;
timestamp = timestamp - timeshift; % Shift to main data set

extension = '12-08-2013';
% extension = '12-25-2013';

sites_all = [2:8 10 20 19 17:-1:11 21:29 40:-1:31 41:50 60:-1:56 54 53 52 51];
sites_colored = [22:24 21 2:4]; % BTC vs. IGF
% sites_colored = [22:24 21]; % only IGF
% sites_colored = [2:4 21]; % only BTC

pcs = [2 3];

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 14;

setFigure(f1,xfac,yfac,fontsize)

hold on

c_signal = [];
celltype = [];

for isite = sites_all
    [scores tmpsignal] = fPCA(isite,extension,timeshift);
    c_signal = [c_signal tmpsignal];
    celltype = [celltype ones(1,size(tmpsignal,2))*isite];
    plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
end

ncolor = 201;
colmap = jet(ncolor);
dose_range = [0 100];
dose_inds = 10.^linspace(log10(max([dose_range(1) .1])),log10(dose_range(2)),floor(ncolor/2));

firstprop = siteprop(sites_colored(1),extension);
for isite = sites_colored
    scores = fPCA(isite,extension,timeshift);
    sprop = siteprop(isite,extension);
    
    colfac = 2*(sprop.lig_index == firstprop.lig_index)-1;
    [tmp colind] = min(abs(sprop.lig_dose - dose_inds));
    mycolor = colmap(ceil(ncolor/2) + colfac*colind,:);
    
    plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',mycolor)
%     plotEllipsis(scores(pcs(1),:),scores(pcs(2),:),mycolor,.5);
    
end

ylabel(['PC ' num2str(pcs(2))])
xlabel(['PC ' num2str(pcs(1))])

clim = log10(dose_range);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(colmap)
colorbar('YTick',log10([1 10 100]),'YTickLabel',{sprop.lig_name,'No Stim',firstprop.lig_name}) % Vertical colorbar

% return
% close all

f2 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f2,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 6;
ncols = 10;

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite))
    plot(repmat(timestamp,1,sum(celltype == isite)),c_signal(:,celltype == isite),'Color',[.7 .7 .7])
    hold on
    plot(timestamp,nanmean(c_signal(:,celltype == isite),2),'color','k','LineWidth',2)
    set(gca,'XLim',[50 200])
    set(gca,'YLim',[-.01 .01])
    plot([120 120],get(gca,'YLim'),'b--')
    s = siteprop(isite,extension);
    titstr = s.lig_name;
    titstr = sprintf('%s %i',titstr,s.lig_dose);
    titstr = sprintf('%s - %s',titstr,s.celltype);
    title(titstr)
end