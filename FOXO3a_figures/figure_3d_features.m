% Figure 3c: pulsatile late response in 184A1
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

% puls_thres = .1;
puls_thres = 0.3;

load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC')

ylabels = {'','','log_{10} Cyt/Nuc','min'};

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

figure

myfeat = [2 3 4 6];

rowstocols = 0.5;
nrows = ceil(length(myfeat)^rowstocols);
ncols = ceil(length(myfeat) / nrows);

colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
resort2 = [6 2 3 4 5 1];
myslopes_feat = nan(length(highdoses),max(myfeat));
mystd_feat = nan(length(highdoses),max(myfeat));
for iplot = 1:length(myfeat)
    icol = myfeat(iplot);

    subplot(nrows,ncols,iplot)
    hold on
    mycolor = lines(nrows);
    legh = [];
    legstr = cell(length(highdoses),1);
    for irow = 1:length(highdoses)
        isite = highdoses(resort2(irow));
        sprop = siteprop(isite);
        legstr{isite == highdoses} = sprop.lig_name(1:3);
        myind = ~isnan(medians(:,irow,icol));
        legh = [legh plot(0:sum(myind)-1,medians(myind,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)];

    %     [axb s] = polyfitZero(1:sum(myind)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    %     legh = [legh plot(0:sum(myind)-1,(0:sum(myind)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)];

        [axb s] = polyfit(0:sum(myind)-1,medians(myind,irow,icol)',1);
        plot(0:sum(myind)-1,(0:sum(myind)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);

        plot([0 sum(myind)-1],[nanmean(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1) nanmean(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1)],'k--','LineWidth',2)
        myslopes_feat(irow,icol) = axb(1);
        Rinv = inv(s.R);
        covmat = sqrt((Rinv*Rinv')*s.normr^2/(s.df));
        mystd_feat(irow,icol) = covmat(1);
    end
    title(features{icol} )
    set(gca,'XLim',[-.5 length(possible_doses)-.5])
    set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
    xlabel('Ligand dose [ng/ml]')
    ylabel(ylabels{iplot})
    
end

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

myalpha_feat = atan(myslopes_feat);
mystdalpha_feat = 1./(1+myslopes_feat.^2) .* mystd_feat;

% myttest_feat = sqrt(s.df)*myslopes_feat./mystd_feat;
% mytestthres_feat = icdf('t',.95,s.df-1);

figure

rankorderqty = myslopes_feat;
rankorderstd = mystd_feat;
titstr = 'Slope';
% rankorderqty = myalpha_feat;
% rankorderstd = mystdalpha_feat;
% titstr = 'Angle';

for iplot = 1:length(myfeat)
    icol = myfeat(iplot);

    subplot(nrows,ncols,iplot)
    hold on

    legh = [];
    for i = 1:length(highdoses)
        legh = [legh errorbar(i,rankorderqty(i,icol),2*rankorderstd(i,icol),markers{i},'Color',colmap(i,:),'LineWidth',2)];
    end
    set(gca,'XTick',1:length(highdoses),'XTickLabel',legstr)
    title(features{icol} )
    ylabel(titstr)
    plot(get(gca,'XLim'),[0 0],'k--','LineWidth',2)
end

h = legend(legh,legstr);

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
myslopes = nan(1,length(highdoses));
mystd = nan(1,length(highdoses));
for irow = 1:length(highdoses)
    sprop = siteprop(highdoses(irow));
    legstr{irow} = sprop.lig_name(1:3);
    legh = [legh plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)];

    [axb s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)
    
    myslopes(irow) = axb(1);
    Rinv = inv(s.R);
    covmat = sqrt((Rinv*Rinv')*s.normr^2/(s.df+1));
    mystd(irow) = covmat(1);
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

myalpha = atan(myslopes);
mystdalpha = 1./(1+myslopes.^2) .* mystd;

% myttest = sqrt((s.df+1))*myslopes./mystd;
% mytestthres = icdf('t',.95,s.df);

figure
hold on

rankorderqty = myslopes;
rankorderstd = mystd;
titstr = 'Slope fraction cells';
% rankorderqty = myalpha;
% rankorderstd = mystdalpha;
% titstr = 'Angle fraction cells';

legh = [];
for i = 1:length(highdoses)
    legh = [legh errorbar(i,rankorderqty(i),2*rankorderstd(i),markers{i},'Color',colmap(i,:),'LineWidth',2)];
end
set(gca,'XTick',1:length(highdoses),'XTickLabel',legstr)
title(titstr)
plot(get(gca,'XLim'),[0 0],'k--','LineWidth',2)

h = legend(legh,legstr,'Location','NorthWest');