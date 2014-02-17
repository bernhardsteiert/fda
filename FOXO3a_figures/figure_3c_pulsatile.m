% Figure 3c: pulsatile late response in 184A1
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70]; % Without FGF
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

puls_thres = .5;

load('./Workspaces/scores_puls')

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,length(features));
highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
    medians(doseind,resort(sprop.lig_index),1:7) = nanmedian(scores_puls(celltypes == isite,:),1);
    medians(doseind,resort(sprop.lig_index),8) = sum(scores_puls(celltypes == isite,1) > puls_thres)/sum(celltypes == isite);
end

resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

nrows = 1;
ncols = 3;

figure

subplot(nrows,ncols,1)

hold on
icol = 8;
legh = [];
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    legh = [legh plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)];
end
title('Fraction cells [%]')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
set(gca,'YTick',0:.05:.40,'YTickLabel',0:5:40)
set(gca,'YLim',[0 .40])
xlabel('Ligand dose [ng/ml]')

subplot(nrows,ncols,2)
hold on
icol = 4;
mycolor = lines(nrows);
legh = [];
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    legh = [legh plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)];
end
title([features{icol} ' [log_{10} Cyt/Nuc]'])
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')


s5 = subplot(nrows,ncols,3);
hold on
icol = 6;
mycolor = lines(nrows);
legh = [];
legstr = cell(length(highdoses),1);
resort2 = [6 2 3 4 5 1];
for irow = 1:length(highdoses)
    isite = highdoses(resort2(irow));
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    legh(irow) = plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6);

    [axb s] = polyfit(1:length(possible_doses)-1,medians(2:end,irow,icol)',1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);
end

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

title([features{icol} ' [min]'])
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')

% subplotpos = get(s5,'Position');
% set(gca,'CLim',[0 1])
% colormap(colmap)
% colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0],'Position',[subplotpos(1)+subplotpos(3) subplotpos(2) .01 subplotpos(4)],'units','normalized') % Vertical colorbar




figure

subplot(nrows,ncols,1)
hold on

% Schematic of how fraction is determined
nnorm = 201;
xnorm = linspace(-5,5,nnorm);
ynorm = normpdf(xnorm,-1.5,1)+.5*normpdf(xnorm,2,1);
ycut = .5; % at .5 when looked at plot(xnorm,ynorm)
% With distribution shifted up
% plot(linspace(inlay_x1,inlay_x2,nnorm),ylim(1) + range(ylim)*(inlay_yscale1+.035) + ynorm/max(ynorm)*range(ylim)*(inlay_yscale2-inlay_yscale1-.05),'k')
% plot(inlay_x1+(ycut-min(xnorm))/range(xnorm)*(inlay_x2-inlay_x1)*[1 1],ylim(1) + range(ylim)*[inlay_yscale1+.035 inlay_yscale2])
% plot([inlay_x1 inlay_x2],ylim(1) + range(ylim)*(inlay_yscale1+.035) + [0 0])
% With distribution in full box
plot(xnorm,ynorm,'k')
ylim = [0 .45];
plot([ycut ycut],ylim)

xlabel('pulsatory strength')
ylabel('distribution')

set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(nrows,ncols,2)
hold on

% Schematic of how amplitude is determined
xpuls = linspace(.5,3.8*pi,nnorm);
ypuls = -cos(xpuls) + xpuls*0.05 - sin(xpuls*.25);
ypuls = ypuls-min(ypuls)+.2*range(ypuls);
[miny minindx] = min(ypuls);
[maxy maxindx] = max(ypuls);
plot(xpuls,ypuls,'k')
plot(xpuls([minindx maxindx]),[miny miny])
plot(xpuls([minindx maxindx]),[maxy maxy])

set(gca,'XLim',[xpuls(1) xpuls(end)])

axPos = get(gca,'Position');
xMinMax = get(gca,'XLim');
yMinMax = get(gca,'YLim');
xAnnotation = axPos(1) + ((mean(xpuls([minindx maxindx]))*[1 1] - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
yAnnotation = axPos(2) + (([miny maxy] - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
annotation('doublearrow',xAnnotation,yAnnotation,'Color','b')

xlabel('time')
ylabel('signal')

set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(nrows,ncols,3)
hold on

% Schematic of how frequency is determined
plot(xpuls,ypuls,'k')
[tmp minindx] = min(abs(xpuls-3));
miny = ypuls(minindx);
plot(xpuls([minindx minindx]),[miny maxy])
plot(xpuls([maxindx maxindx]),[miny maxy])

set(gca,'XLim',[xpuls(1) xpuls(end)])

axPos = get(gca,'Position');
xMinMax = get(gca,'XLim');
yMinMax = get(gca,'YLim');
xAnnotation = axPos(1) + ((xpuls([minindx maxindx]) - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
yAnnotation = axPos(2) + ((mean([miny maxy])*[1 1] - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
annotation('doublearrow',xAnnotation,yAnnotation,'Color','b')

xlabel('time')
ylabel('signal')

set(gca,'XTick',[])
set(gca,'YTick',[])