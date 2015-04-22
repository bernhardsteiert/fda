% Figure 2a: Linear correlation of pERK/pAKT vs. pS294/pS253
addpath('./Functions/')

[westerndata,description] = xlsread('./Workspaces/westernBlotData');
liglabels = {'NS','EGF','IGF','FGF','HRG','HGF','EPR','BTC'};
exind = 1;
timeind = 2;
ligind = 3;
obsind = 4:10;
myobs = [4 5 8 9];
normind = 11;
uni_ligs = unique(westerndata(:,ligind));
resort = [1 3 5 6 2 8 7 4];
uni_ligs = uni_ligs(resort);
liglabels = liglabels(resort);

experiments_investigated = 2:3; % Write 2:3 to ignore experiment 1

westerndata = westerndata(ismember(westerndata(:,exind),experiments_investigated),:);
data = westerndata;
data(:,obsind) = westerndata(:,obsind)./repmat(westerndata(:,normind),1,length(obsind));
for iexp = experiments_investigated
    tmpdata = data(data(:,exind) == iexp,obsind);
    data(data(:,exind) == iexp,obsind) = data(data(:,exind) == iexp,obsind)./repmat(max(tmpdata,[],1),sum(data(:,exind) == iexp),1);
end
data(:,obsind) = log2(data(:,obsind));

rowstocols = 0.5;
nrows = ceil(length(obsind)^rowstocols);
ncols = ceil(length(obsind) / nrows);

% Plotting raw data
% close all
% colmap = lines(length(experiments_investigated));
% for ifig = 1:length(uni_ligs)
%     figure
%     for iplot = 1:length(obsind)
%         subplot(nrows,ncols,iplot)
%         hold on
%         for iexp = experiments_investigated
%             myind = data(:,exind) == iexp & data(:,ligind) == uni_ligs(ifig);
%             plot(data(myind,timeind),data(myind,obsind(iplot)),'o-','Color',colmap(iexp,:))
%             title([liglabels{ifig} ' - ' description{obsind(iplot)}])
%         end
%     end
% end

% Scaling of experiments
scaledLigs = data(data(:,exind) == experiments_investigated(1),ligind);
scaledTime = data(data(:,exind) == experiments_investigated(1),timeind);
pinit = data(data(:,exind) == experiments_investigated(1),myobs);
pinit = [ones(length(experiments_investigated)-1,length(myobs)); pinit];
[optRes,~,~,~,~,~,jacobian] = lsqnonlin(@(p) res(p,data(:,[exind timeind myobs])),pinit);

Cov=inv(jacobian'*jacobian);
scaledData = optRes(length(experiments_investigated):end,:);
scaledStd = sqrt(diag(Cov));
scaledStd = scaledStd((length(experiments_investigated)-1)*length(myobs)+1:end);
scaledStd = reshape(scaledStd,size(scaledData));

rowstocols = 0.5;
nrows = ceil(size(scaledData,2)^rowstocols);
ncols = ceil(size(scaledData,2) / nrows);

% Plotting time-courses
colmap = [linspace(0,1,length(uni_ligs)-1)' ones(length(uni_ligs)-1,1) ones(length(uni_ligs)-1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
colmap = [colmap; [0 0 0]];
% close all
% figure
% for iplot = 1:size(scaledData,2)
%     subplot(nrows,ncols,iplot)
%     hold on
%     colcount = 1;
%     for ilig = uni_ligs(2:end)'
%         myind = logical([1; scaledLigs(2:end) == ilig]); % Plot common timepoint zero for all ligands
%         if sum(~isnan(scaledData(myind,iplot))) > 1 % not only NaN
%             plot(scaledTime(myind),scaledData(myind,iplot),'o-','Color',colmap(colcount,:))
%             errorbar(scaledTime(myind),scaledData(myind,iplot),scaledStd(myind,iplot)./2,'Color',colmap(colcount,:))
%             colcount = colcount + 1;
%         end
%     end
%     title(description{myobs(iplot)})
%     ylabel('log_2 fold change [au]')
%     xlabel('time [min]')
%     set(gca,'XLim',[-10 490])
% end

figure
legh = [];
for iplot = 1:size(scaledData,2)
    subplot(nrows,ncols,iplot)
    hold on
    colcount = 1;
    for ilig = uni_ligs(2:end)'
        myind = logical([1; scaledLigs(2:end) == ilig]); % Plot common timepoint zero for all ligands
        if sum(~isnan(scaledData(myind,iplot))) > 1 % not only NaN
            tmpx = [scaledTime(myind); flipud(scaledTime(myind))];
            tmpy = [scaledData(myind,iplot) + scaledStd(myind,iplot)./2; flipud(scaledData(myind,iplot) - scaledStd(myind,iplot)./2)];
            
            ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
            set(ltmp, 'FaceColor', colmap(colcount,:)*0.1+0.9, 'EdgeColor', colmap(colcount,:)*0.1+0.9);
            ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
            set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(colcount,:)*0.3+0.7);
                                    
                                    
            legh = [legh plot(scaledTime(myind),scaledData(myind,iplot),'x-','Color',colmap(colcount,:))];
%             errorbar(scaledTime(myind),scaledData(myind,iplot),scaledStd(myind,iplot)./2,'Color',colmap(colcount,:))
            colcount = colcount + 1;
        end
    end
    title(description{myobs(iplot)})
    ylabel('log_2 fold change [au]')
    xlabel('time [min]')
    set(gca,'XLim',[-10 490])
     set(gca,'YLim',[-7 1],'YTick',[-6:2:0])
end

legend(legh,liglabels{2:end-1})

% --------

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

puls_thres = 0.3;

late_ws = load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC');
early_ws = load('./Workspaces/scores_early_5basis_noFGF_newBTC');

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,1);
highdoses = [];

rng(0) % Make sure that bootstrap samples are reproducible
% ratio_fun = @(x) sum(x > puls_thres) / length(x);

for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    medians(doseind,resort(sprop.lig_index),1:5) = nanmedian(early_ws.scores_early(:,early_ws.celltypes == isite),2);
    medians(doseind,resort(sprop.lig_index),6) = sum(late_ws.scores_puls(late_ws.celltypes == isite,1) > puls_thres)/sum(late_ws.celltypes == isite);
end
medians(1,5,:) = nanmean(medians(1,:,:),2);


close all
% Plotting correlation diagram
markers = {'o','s','v','d','^','>','<'};
legh = [];
figure
colcount = 1;
erkaktratio = [];
erkaktratio_std = [];
foxositesratio = [];
foxositesratio_std = [];

for ilig = uni_ligs(2:end)' % Ignore NS case
    myind = scaledLigs == ilig;
    if sum(~isnan(scaledData(myind,:))) % not only NaN
        erkaktratio = [erkaktratio scaledData(myind,1)-scaledData(myind,2)];
        foxositesratio = [foxositesratio scaledData(myind,3)-scaledData(myind,4)];
        erkaktratio_std = [erkaktratio_std sqrt(sum(scaledStd(myind,3:4).^2,2))];
        foxositesratio_std = [foxositesratio_std sqrt(sum(scaledStd(myind,1:2).^2,2))];
%         plot(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4)];
%         legh = [legh errorbar(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4),sqrt(sum(scaledStd(myind,3:4).^2./2,2)),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:))];
%         h = herrorbar(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4),sqrt(sum(scaledStd(myind,1:2).^2./2,2)));
        subplot(3,4,1)
        hold on
        legh = [legh plot(mean(scaledData(myind,1)-scaledData(myind,2)),mean(scaledData(myind,3)-scaledData(myind,4)),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12)];
        xlabel(['log_{2} ' description{myobs(1)} '/' description{myobs(2)}])
        ylabel(['log_{2} ' description{myobs(3)} '/' description{myobs(4)}])
        subplot(3,4,2)
        hold on
        plot(mean(scaledData(myind,1)),mean(scaledData(myind,3)),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12);
        xlabel(['log_{2} ' description{myobs(1)}])
        ylabel(['log_{2} ' description{myobs(3)}])
        subplot(3,4,3)
        hold on
        plot(mean(scaledData(myind,2)),mean(scaledData(myind,4)),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12);
        xlabel(['log_{2} ' description{myobs(2)}])
        ylabel(['log_{2} ' description{myobs(4)}])
        for i = 1:3
            subplot(3,4,3+i)
            hold on
            myind2 = find(myind);
            plot(scaledData(myind2(1),1),medians(end,size(erkaktratio,2),i),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12);
            xlabel(sprintf('log_{2} %s (t = %d)',description{myobs(1)},scaledTime(myind2(1))))
            ylabel(sprintf('Early PC%d',i))
        end
        subplot(3,4,7)
        hold on
        myind2 = find(myind);
        plot(scaledData(myind2(1),1)-scaledData(myind2(1),2),medians(end,size(erkaktratio,2),6),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12);
        xlabel(sprintf('log_{2} %s/%s (t = %d)',description{myobs(1)},description{myobs(2)},scaledTime(myind2(1))))
        ylabel('Fraction of pulsing cells')
        for i = 1:3
            subplot(3,4,7+i)
            hold on
            plot(mean(scaledData(myind,2)),medians(end,size(erkaktratio,2),i),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12);
            xlabel(['log_{2} ' description{myobs(2)}])
            ylabel(sprintf('Early PC%d',i))
        end
        subplot(3,4,11)
        hold on
        plot(mean(scaledData(myind,2)),medians(end,size(erkaktratio,2),6),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:),'MarkerSize',12);
        xlabel(['log_{2} ' description{myobs(2)}])
        ylabel('Fraction of pulsing cells')
%         for i = 1:3
%             figure(6+i)
%             plot(mean(scaledData(myind,2)),medians(end,size(erkaktratio,2),i),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:));
%             xlabel(['log_{2} ' description{myobs(2)}])
%             ylabel(sprintf('Early PC%d',i))
%         end

%         set(h,'Color',colmap(colcount,:))
%         set(h(2),'LineStyle', 'none')
%         errorbar_tick(legh(end),0)
        colcount = colcount + 1;
    end
end
subplot(3,4,1)
axb = lsqnonlin(@(axb) [(erkaktratio(:)-(foxositesratio(:)-axb(2))./axb(1))./(erkaktratio_std(:)); (foxositesratio(:)-axb(1)*foxositesratio(:)-axb(2))./(foxositesratio_std(:))],[1 0]);
plot([min(min(erkaktratio)) max(max(erkaktratio))],[min(min(erkaktratio)) max(max(erkaktratio))]*axb(1) + axb(2),'k--','LineWidth',2)
chi2 = sum(([(erkaktratio(:)-(foxositesratio(:)-axb(2))./axb(1))./(erkaktratio_std(:)); (foxositesratio(:)-axb(1)*foxositesratio(:)-axb(2))./(foxositesratio_std(:))]).^2);
% 
% legend(legh,liglabels{2:end},'Location','NorthWest')
% xlabel(['log_{2} ' description{myobs(1)} '/' description{myobs(2)}])
% ylabel(['log_{2} ' description{myobs(3)} '/' description{myobs(4)}])
% 
set(gca,'XLim',[-7.5 5],'YLim',[-3 3.5],'XTick',[-6:2:4])

% text(1,-1,['\chi^2/N = ' num2str(chi2/(2*length(erkaktratio(:))))])

