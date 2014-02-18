% Figure 3c: pulsatile late response in 184A1
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70]; % Without FGF
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

puls_thres = .5;

load('./Workspaces/scores_early')

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,3);
highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
    medians(doseind,resort(sprop.lig_index),:) = nanmedian(scores_early(:,celltypes == isite),2);
end

resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

nrows = 1;
ncols = 3;

figure

setFigure(gcf,1,.5,12)

subplot(nrows,ncols,1)

hold on
icol = 1;
legh = [];
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    [axb s] = polyfit(0:length(possible_doses)-1,medians(1:end,irow,icol)',1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);
end
title('Harmonic 1')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% set(gca,'YTick',0:.05:.40,'YTickLabel',0:5:40)
% set(gca,'YLim',[0 .40])
xlabel('Ligand dose [ng/ml]')


subplot(nrows,ncols,2)
hold on
icol = 2;
mycolor = lines(nrows);
legh = [];
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    legh = [legh plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)];
end
title('Harmonic 2')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')


s5 = subplot(nrows,ncols,3);
hold on
icol = 3;
mycolor = lines(nrows);
legh = [];
legstr = cell(length(highdoses),1);
resort2 = [6 2 3 4 5 1];
for irow = 1:length(highdoses)
    isite = highdoses(resort2(irow));
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    legh(irow) = plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6);

    [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2);
end

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

title('Harmonic 3')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')

figure

setFigure(gcf,1,.5,12)

subplot(nrows,ncols,1)

hold on
icol = 1;
legh = [];
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    % Fill all parameters
    qFit = logical([1 1 1 1]);
    pFix = [];
    pinit = [medians(1,irow,icol) medians(end,irow,icol) length(possible_doses)/2 .01];
    optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf -Inf -Inf],[Inf Inf Inf Inf],optimset('Display','off'));
    
    % Plotting
    p = nan(size(qFit));
    p(qFit) = optRes;
    p(~qFit) = pFix;
    legh = [legh plot(linspace(0,length(possible_doses)-1,201),f(p,linspace(0,length(possible_doses)-1,201)),'-','Color',colmap(irow,:),'LineWidth',2)];
end
title('Harmonic 1')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% set(gca,'YTick',0:.05:.40,'YTickLabel',0:5:40)
% set(gca,'YLim',[0 .40])
xlabel('Ligand dose [ng/ml]')


subplot(nrows,ncols,2)
hold on
icol = 2;
mycolor = lines(nrows);
legh = [];
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    % Fill all parameters
    qFit = logical([1 1 1 1]);
    pFix = [];
    pinit = [medians(1,irow,icol) medians(end,irow,icol) length(possible_doses)/2 .01];
    optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf -Inf -Inf],[Inf Inf Inf Inf],optimset('Display','off'));
    
    % Plotting
    p = nan(size(qFit));
    p(qFit) = optRes;
    p(~qFit) = pFix;
    legh = [legh plot(linspace(0,length(possible_doses)-1,201),f(p,linspace(0,length(possible_doses)-1,201)),'-','Color',colmap(irow,:),'LineWidth',2)];
end
title('Harmonic 2')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')


s5 = subplot(nrows,ncols,3);
hold on
icol = 3;
mycolor = lines(nrows);
legh = [];
legstr = cell(length(highdoses),1);
resort2 = [6 2 3 4 5 1];
for irow = 1:length(highdoses)
    isite = highdoses(resort2(irow));
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    legh(irow) = plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6);

    % Fill all parameters
    qFit = logical([1 1 1 1]);
    pFix = [];
    pinit = [medians(1,irow,icol) medians(end,irow,icol) length(possible_doses)/2 .01];
    optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf -Inf -Inf],[Inf Inf Inf Inf],optimset('Display','off'));
    
    % Plotting
    p = nan(size(qFit));
    p(qFit) = optRes;
    p(~qFit) = pFix;
    plot(linspace(0,length(possible_doses)-1,201),f(p,linspace(0,length(possible_doses)-1,201)),'-','Color',colmap(irow,:),'LineWidth',2);
end

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

title('Harmonic 3')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')

sites_colored = [11:17 64:70];

pcs = [2 3];

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 12;

setFigure(f1,xfac,yfac,fontsize)

plot(scores_early(2,~ismember(celltypes,sites_colored)),scores_early(3,~ismember(celltypes,sites_colored)),'o','MarkerSize',1,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none')
hold on

ncolor = 201;
colmap = jet(ncolor);
dose_range = [0 100];
dose_inds = 10.^linspace(log10(max([dose_range(1) 1])),log10(dose_range(2)),floor(ncolor/2));

firstprop = siteprop(sites_colored(1));
for isite = sites_colored
    sprop = siteprop(isite);
    
    colfac = 2*(sprop.lig_index == firstprop.lig_index)-1;
    [tmp colind] = min(abs(sprop.lig_dose - dose_inds));
    mycolor = colmap(ceil(ncolor/2) + colfac*colind,:);
    
    plot(scores_early(pcs(1),celltypes == isite),scores_early(pcs(2),celltypes == isite),'.','Color',mycolor)
%     plotEllipsis(scores(pcs(1),:),scores(pcs(2),:),mycolor,.5);
    
end

xlim = [-.15 .25];
ylim = [-.08 .15];
set(gca,'XLim',xlim,'YLim',ylim)

% set(gca,'XLim',xlim)

% axisEqual(get(gcf,'Position'))

ylabel('score of transient harmonic')
xlabel('score of sustained harmonic')

% ylabel(['PC ' num2str(pcs(2))])
arrow([-.13 -.06],[.2 -.06],'Width',.5,'Length',7)
% set(gca,'YTick',-.1:.1:.2)
% xlabel(['PC ' num2str(pcs(1))])
arrow([-.13 .05],[.02 .14],'Width',.5,'Length',7)

clim = log10(dose_range);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(colmap)
colorbar('YTick',log10([1 10 100]),'YTickLabel',{'BTC','No Stim','IGF'}) % Vertical colorbar

