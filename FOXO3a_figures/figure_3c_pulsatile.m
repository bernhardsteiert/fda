% Figure 3c: pulsatile late response in 184A1
close all;
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

% puls_thres = .1;
puls_thres = 0.3;

% load('./Workspaces/scores_puls')
load('./Workspaces/scores_puls_corrected_retracked_all_cleaned')

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,length(features)+1);
highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
%     medians(doseind,resort(sprop.lig_index),1:7) = nanmedian(scores_puls(celltypes == isite,:),1);
    medians(doseind,resort(sprop.lig_index),1:7) = nanmean(scores_puls(celltypes == isite & scores_puls(:,1)' > puls_thres,:),1);
    medians(doseind,resort(sprop.lig_index),8) = sum(scores_puls(celltypes == isite,1) > puls_thres)/sum(celltypes == isite);
end
medians(1,5,:) = nanmean(medians(1,:,:),2);

resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

nrows = 1;
ncols = 3;

figure

hold on
icol = 8;
legh = [];
legstr = cell(length(highdoses),1);
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
for irow = 1:length(highdoses)
    sprop = siteprop(highdoses(irow));
    legstr{irow} = sprop.lig_name(1:3);
    legh = [legh plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)];

    [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)
end
title('Fraction cells [%]')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
set(gca,'YTick',0:.1:1,'YTickLabel',0:10:100)
set(gca,'YLim',[0 1])
xlabel('Ligand dose [ng/ml]')

% subplot(nrows,ncols,2)
% hold on
% icol = 2;
% mycolor = lines(nrows);
% legh = [];
% for irow = 1:length(highdoses)
%     plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)
% 
% %     [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
% %     legh = [legh plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)];
% 
%     [axb s] = polyfit(0:length(possible_doses)-1,medians(:,irow,icol)',1);
%     plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);
%     
%     plot([0 length(possible_doses)-1],[nanmean(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1) nanmean(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1)],'k--','LineWidth',2)
% end
% title([features{icol} ' [log_{10} Cyt/Nuc]'])
% set(gca,'XLim',[-.5 length(possible_doses)-.5])
% set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% xlabel('Ligand dose [ng/ml]')
% 
% 
% s5 = subplot(nrows,ncols,3);
% hold on
% icol = 6;
% mycolor = lines(nrows);
% legh = [];
% legstr = cell(length(highdoses),1);
% resort2 = [6 2 3 4 5 1];
% for irow = 1:length(highdoses)
%     isite = highdoses(resort2(irow));
%     sprop = siteprop(isite);
%     legstr{isite == highdoses} = sprop.lig_name(1:3);
%     legh(irow) = plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6);
% 
% %     [axb s] = polyfit(1:length(possible_doses)-1,medians(2:end,irow,icol)',1);
%     [axb s] = polyfit(0:length(possible_doses)-1,medians(:,irow,icol)',1);
%     plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);
%     
%     plot([0 length(possible_doses)-1],[nanmedian(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1) nanmedian(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1)],'k--','LineWidth',2)
% end

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end
% 
% title([features{icol} ' [min]'])
% set(gca,'XLim',[-.5 length(possible_doses)-.5])
% set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% xlabel('Ligand dose [ng/ml]')

% figure
% 
% subplot(nrows,ncols,1)
% 
% hold on
% icol = 8;
% legh = [];
% colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
% colmap = hsv2rgb(colmap(1:end-1,:));
% markers = {'o','s','v','d','^','>'};
% f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));
% for irow = 1:length(highdoses)
%     plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)
%     
%     % Fix lower asymptotic to mean(medians)
% %     qFit = logical([0 1 1 1]);
% %     pFix = mean(medians(1,:,icol));
% %     pinit = [medians(end,irow,icol) length(possible_doses)/2 .1];
% %     optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf .2],[Inf Inf Inf],optimset('Display','off'));
% 
%     % Fill all parameters
%     qFit = logical([1 1 1 1]);
%     pFix = [];
%     pinit = [medians(1,irow,icol) medians(end,irow,icol) length(possible_doses)/2 .1];
%     optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf -Inf .2],[Inf Inf Inf Inf],optimset('Display','off'));
%     
%     % Plotting
%     p = nan(size(qFit));
%     p(qFit) = optRes;
%     p(~qFit) = pFix;
%     legh = [legh plot(linspace(0,length(possible_doses)-1,201),f(p,linspace(0,length(possible_doses)-1,201)),'-','Color',colmap(irow,:),'LineWidth',2)];
% end
% title('Fraction cells [%]')
% set(gca,'XLim',[-.5 length(possible_doses)-.5])
% set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% set(gca,'YTick',0:.05:.40,'YTickLabel',0:5:40)
% set(gca,'YLim',[0 .40])
% xlabel('Ligand dose [ng/ml]')
% 
% subplot(nrows,ncols,2)
% hold on
% icol = 4;
% mycolor = lines(nrows);
% legh = [];
% for irow = 1:length(highdoses)
%     plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)
% 
%     % Fix lower asymptotic to mean(medians)
%     qFit = logical([0 1 1 1]);
%     pFix = mean(medians(1,:,icol));
%     pinit = [medians(end,irow,icol) length(possible_doses)/2 .1];
%     optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf .2],[Inf Inf Inf],optimset('Display','off'));
% 
%     % Fill all parameters
% %     qFit = logical([1 1 1 1]);
% %     pFix = [];
% %     pinit = [medians(1,irow,icol) medians(end,irow,icol) length(possible_doses)/2 .1];
% %     optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf -Inf .2],[Inf Inf Inf Inf],optimset('Display','off'));
%     
%     % Plotting
%     p = nan(size(qFit));
%     p(qFit) = optRes;
%     p(~qFit) = pFix;
%     legh = [legh plot(linspace(0,length(possible_doses)-1,201),f(p,linspace(0,length(possible_doses)-1,201)),'-','Color',colmap(irow,:),'LineWidth',2)];
% end
% title([features{icol} ' [log_{10} Cyt/Nuc]'])
% set(gca,'XLim',[-.5 length(possible_doses)-.5])
% set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% xlabel('Ligand dose [ng/ml]')
% 
% 
% s5 = subplot(nrows,ncols,3);
% hold on
% icol = 6;
% mycolor = lines(nrows);
% legh = [];
% legstr = cell(length(highdoses),1);
% resort2 = [6 2 3 4 5 1];
% for irow = 1:length(highdoses)
%     isite = highdoses(resort2(irow));
%     sprop = siteprop(isite);
%     legstr{isite == highdoses} = sprop.lig_name(1:3);
%     legh(irow) = plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6);
% 
%     [axb s] = polyfit(1:length(possible_doses)-1,medians(2:end,irow,icol)',1);
%     plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);
% end
% 
% h = legend(legh,legstr);
% ch = get(h,'child');
% for ileg = 1:length(ch)/3
%     ilegch = (ileg-1)*3+2;
%     set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
% end
% 
% title([features{icol} ' [min]'])
% set(gca,'XLim',[-.5 length(possible_doses)-.5])
% set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
% xlabel('Ligand dose [ng/ml]')


figure

subplot(nrows,ncols,1)
hold on

% Schematic of how fraction is determined
nnorm = 201;
xnorm = linspace(-5,5,nnorm);
ynorm = normpdf(xnorm,-1.5,1)+.5*normpdf(xnorm,2,1);
ycut = .5; % at .5 when looked at plot(xnorm,ynorm)
plot(xnorm,ynorm,'k')
ylim = [0 .45];
plot([ycut ycut],ylim)

xlabel('pulsatory strength')
ylabel('probability density')

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


figure
hold on

kswidth = 20;

unstim = [10 11 31 50 51];
for i = 2:length(unstim)
    siteprop(unstim(i))
    celltypes(celltypes == unstim(i)) = unstim(1);
end
highdoses = [highdoses unstim(1)];
[dists_boxcox,lambda,c] = boxcox(scores_puls(:,1));
thres_trafo = boxcox_apply(puls_thres,lambda,c);

legh = [];
legstr = {};
colmap = [colmap; [0 0 0]];
for isite = highdoses([1 6 5 7])
    s = siteprop(isite);

    [f,xi] = ksdensity(dists_boxcox(celltypes == isite),'width',range(dists_boxcox)./kswidth);
    legh = [legh plot(xi,f,'Color',colmap(isite == highdoses,:))];

    legstr{end+1} = s.lig_name;
    
end
legstr{end} = 'NS';
set(gca,'XLim',[min(dists_boxcox) max(dists_boxcox)]+[-.2 .2]*range(dists_boxcox))
plot([thres_trafo thres_trafo],get(gca,'YLim'),'k--')
% set(gca,'XLim',[-6 2])
% title(s.celltype)
legend(legh,legstr)

xlabel('pulsatory score')
ylabel('probability density')