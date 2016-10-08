% Figure 3c: pulsatile late response in 184A1
close all;
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

% puls_thres = .1;
puls_thres = 0.3;

% load('./Workspaces/scores_puls')
load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC')
late_ws = load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC');
early_ws = load('./Workspaces/scores_early_5basis_noFGF_newBTC');

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,length(features)+1+2);
highdoses = [];

rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x) sum(x > puls_thres) / length(x);

for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
%     medians(doseind,resort(sprop.lig_index),1:7) = nanmedian(scores_puls(celltypes == isite,:),1);
    medians(doseind,resort(sprop.lig_index),1:7) = nanmean(scores_puls(celltypes == isite & scores_puls(:,1)' > puls_thres,:),1);
    medians(doseind,resort(sprop.lig_index),8) = sum(scores_puls(celltypes == isite,1) > puls_thres)/sum(celltypes == isite);
    medians(doseind,resort(sprop.lig_index),9:10) = bootci(2000,{ratio_fun,scores_puls(celltypes == isite,1)},'alpha',.32);
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
    h=errorbar(0:length(possible_doses)-1,medians(:,irow,icol),medians(:,irow,icol)-medians(:,irow,icol+1),medians(:,irow,icol+2)-medians(:,irow,icol),'Color',colmap(irow,:));
    set(h,'linestyle','none')
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

figure

% subplot(nrows,ncols,1)

hold on
icol = 8;
legh = [];
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));
mymean = mean(medians(1,:,8));
for irow = 1:length(highdoses)
    legh = [legh plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor','none','MarkerSize',6)];
    medians(1,irow,icol) = mymean;
    % Fix lower asymptotic to mean(medians)
%     qFit = logical([0 1 1 1]);
%     pFix = mean(medians(1,:,icol));
%     pinit = [medians(end,irow,icol) length(possible_doses)/2 .1];
%     optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf .2],[Inf Inf Inf],optimset('Display','off'));

    % Fill all parameters
    qFit = logical([1 1 1 1]);
    pFix = [];
    pinit = [medians(1,irow,icol) medians(end,irow,icol) length(possible_doses)/2 .1];
    optRes = lsqnonlin(@(p) objFunHill(p,0:length(possible_doses)-1,medians(:,irow,icol)',qFit,pFix),pinit,[-Inf -Inf -Inf .2],[Inf Inf Inf Inf],optimset('Display','off'));
    
    % Plotting
    p = nan(size(qFit));
    p(qFit) = optRes;
    p(~qFit) = pFix;
    plot(linspace(0,length(possible_doses)-1,201),f(p,linspace(0,length(possible_doses)-1,201)),'-','Color',colmap(irow,:),'LineWidth',2)
    
    h=errorbar(0:length(possible_doses)-1,medians(:,irow,icol),medians(:,irow,icol)-medians(:,irow,icol+1),medians(:,irow,icol+2)-medians(:,irow,icol),'Color',colmap(irow,:));
    set(h,'linestyle','none')
end
title('Fraction cells [%]')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
set(gca,'YTick',0:.1:1,'YTickLabel',0:10:100)
set(gca,'YLim',[0 1])
xlabel('Ligand dose [ng/ml]')

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end
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
    siteprop(unstim(i));
    celltypes(celltypes == unstim(i)) = unstim(1);
end
highdoses = [highdoses unstim(1)];
c = 0;
lambda = 0.0500878494830929582581902081984;
dists_boxcox = boxcox_apply(scores_puls(:,1),lambda,c);
thres_trafo = boxcox_apply(puls_thres,lambda,c);
mymin = 0.00360508363121257550953924209125;
mymax = 2.34695221468680559340214131225;
mymin_trafo = boxcox_apply(mymin,lambda,c);
mymax_trafo = boxcox_apply(mymax,lambda,c);
myrange = mymax_trafo - mymin_trafo;

dists_boxcox = (dists_boxcox-mymin_trafo)./myrange;
thres_trafo = (thres_trafo-mymin_trafo)./myrange;

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

figure
hold on

% ------------------------------------------
% Early and late scores are not aligned yet; Sorting is done in the following

sites_all_early = [4:17 31:37 44:57 64:69];
sites_all_late = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69];

sites_all_sorted = sort(sites_all_early);
dists_sorted = nan*dists_boxcox;
early_sorted = nan*early_ws.scores_early;
celltypes_sorted = nan*celltypes;

myind = 1;
for isite = 1:length(sites_all_sorted)
    myind2 = late_ws.celltypes == sites_all_sorted(isite);
    dists_sorted(myind:myind+sum(myind2)-1) = dists_boxcox(myind2);
    early_sorted(:,myind:myind+sum(myind2)-1) = early_ws.scores_early(:,early_ws.celltypes == sites_all_sorted(isite));
    celltypes_sorted(myind:myind+sum(myind2)-1) = late_ws.celltypes(myind2);
    myind = myind+sum(myind2);
end

sites_labels = {'EGF','IGF','HRG','HGF','EPR','BTC'};

highdoses = 1:6;
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
ligand_dose = [100 50 20 10 5 2.5 0];

resort = [2 3 4 1 6 5];
for ilig2 = 1:6
    subplot(3,2,ilig2)
    ilig = resort(ilig2);
    plot(early_sorted(2,:),dists_sorted,'.','Color',[.7 .7 .7])
    hold on
    
    isBTC = ilig == 6;
    if mod(floor((sites_all_sorted((ilig-1)*7+1)-1)/10),2)
        loopind = sites_all_sorted((ilig-1)*7+1):sites_all_sorted(ilig*(7)-isBTC);
    else
        loopind = sites_all_sorted(ilig*(7)-isBTC):-1:sites_all_sorted((ilig-1)*7+1);
    end
    for isite = loopind
        s = siteprop(isite);
        mycol = rgb2hsv(colmap(ilig2,:));
        mycol(3) = find(s.lig_dose == ligand_dose) / length(ligand_dose);
        mycol = hsv2rgb(mycol);
        plot(early_sorted(2,celltypes_sorted == isite),dists_sorted(celltypes_sorted == isite),'.','Color',mycol)
        title(sites_labels{ilig})
    end
    plot(get(gca,'XLim'),[thres_trafo thres_trafo],'k--')
    
    xlabel('Early PC2')
    ylabel('Pulsatory score')
end

%%
figure
hold on
set(gca,'XLim',[-.1 ,.2])
set(gca,'YLim',[0 1],'YTick',0:.1:1,'YTickLabel',0:10:100)
set(gca,'ZLim',[-.01 ,.09])
resort = [4 1 nan 2 3 6 5]; % Relative to platemap
pcs = [2 3];

early_pc_med = nan([size(medians(:,:,1)) size(early_ws.scores_early,1)]);
highdoses = [];
for isite = sites_all_sorted
    s = siteprop(isite);
    if s.lig_dose == 100
        highdoses = [highdoses isite];
    end
    for ipc = 1:size(early_pc_med,3)
        early_pc_med(possible_doses == s.lig_dose,resort(s.lig_index),ipc) = nanmedian(early_ws.scores_early(ipc,early_ws.celltypes == isite));
    end
end
early_pc_med(1,:,:) = repmat(nanmean(early_pc_med(1,:,:),2),1,size(early_pc_med,2),1);
resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

icol = 8;
legh = [];
for irow = 1:length(highdoses)
    sprop = siteprop(highdoses(irow));
    legstr{irow} = sprop.lig_name(1:3);
    legh = [legh plot3(early_pc_med(:,irow,pcs(1)),medians(:,irow,icol),early_pc_med(:,irow,pcs(2)),'--','Color',colmap(irow,:))];
    [err,coeff,p0] = fit_3D_data(early_pc_med(:,irow,pcs(1)),medians(:,irow,icol),early_pc_med(:,irow,pcs(2)),'line','off','off');
    plot3(p0(1)+[-1e3 1e3]*coeff(1),p0(2)+[-1e3 1e3]*coeff(2),p0(3)+[-1e3 1e3]*coeff(3),'Color',colmap(irow,:))
    for idose = 1:size(medians,1)
        mycolor = rgb2hsv(colmap(irow,:));
        mycolor(2:3) = idose/size(medians,1)*mycolor(2:3);
        mycolor = hsv2rgb(mycolor);
        plot3(early_pc_med(idose,irow,pcs(1)),medians(idose,irow,icol),early_pc_med(idose,irow,pcs(2)),markers{irow},'MarkerFaceColor',mycolor,'MarkerEdgeColor',mycolor,'MarkerSize',12);
    end
end
xlabel(sprintf('Early PC%d',pcs(1)))
ylabel('Fraction cells [%]')
zlabel(sprintf('Early PC%d',pcs(2)))
set(gca,'YLim',[0 1],'YTick',0:.1:1,'YTickLabel',0:10:100)

legend(legh,legstr)

%% EGF dose response + AKTi
akti = load('Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC_ATKi');
aktipc2 = load('Workspaces/scores_early_5basis_noFGF_AKTi');
celltypes = [celltypes akti.celltypes];
scores_puls = [scores_puls; akti.scores_puls];
scores_early = [early_ws.scores_early aktipc2.scores_early];
celltypes_early = [early_ws.celltypes akti.celltypes];

sites = [2 4:10];
sites = sites(end:-1:1);
early_pc2 = nan(1,length(sites));
frac_puls = early_pc2;
doses = early_pc2;

figure
hold on
for isite = sites
    s = siteprop(isite);
    doses(isite == sites) = s.lig_dose;
    early_pc2(isite == sites) = nanmedian(scores_early(2,celltypes_early == isite));
    frac_puls(isite == sites) = sum(scores_puls(celltypes == isite,1) > puls_thres)/sum(celltypes == isite);
    
    mycolor = rgb2hsv([0 1 1]);
    mycolor(2:3) = find(isite==sites)/length(sites)*mycolor(2:3);
    mycolor = hsv2rgb(mycolor);
    
    plot(early_pc2(isite == sites),frac_puls(isite == sites),'o','Color',mycolor,'MarkerFaceColor',mycolor)
end
title(s.lig_name)
xlabel('fPC2')
ylabel('Fraction of pulsing cells [%]')

plot(early_pc2(1:end-1),frac_puls(1:end-1),'--','Color',mycolor)