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
data(:,obsind) = log2(westerndata(:,obsind)./repmat(westerndata(:,normind),1,length(obsind)));

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
%             errorbar(scaledTime(myind),scaledData(myind,iplot),scaledStd(myind,iplot),'Color',colmap(colcount,:))
%             colcount = colcount + 1;
%         end
%     end
%     title(description{myobs(iplot)})
%     ylabel('log_2 fold change [au]')
%     xlabel('time [min]')
%     set(gca,'XLim',[-10 490])
% end

% Plotting correlation diagram
markers = {'o','s','p','d','^','o'};
legh = [];
figure
hold on
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
        legh = [legh errorbar(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4),sqrt(sum(scaledStd(myind,3:4).^2,2)),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:))];
        h = herrorbar(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4),sqrt(sum(scaledStd(myind,1:2).^2,2)));
        set(h,'Color',colmap(colcount,:))
        set(h(2),'LineStyle', 'none')
        colcount = colcount + 1;
    end
end
axb = lsqnonlin(@(axb) [(erkaktratio(:)-(foxositesratio(:)-axb(2))./axb(1))./(erkaktratio_std(:)); (foxositesratio(:)-axb(1)*foxositesratio(:)-axb(2))./(foxositesratio_std(:))],[1 0]);
plot([min(min(erkaktratio)) max(max(erkaktratio))],[min(min(erkaktratio)) max(max(erkaktratio))]*axb(1) + axb(2),'k--','LineWidth',2)
chi2 = sum(([(erkaktratio(:)-(foxositesratio(:)-axb(2))./axb(1))./(erkaktratio_std(:)); (foxositesratio(:)-axb(1)*foxositesratio(:)-axb(2))./(foxositesratio_std(:))]).^2);

legend(legh,liglabels{2:end-1},'Location','NorthWest')
xlabel(['log_{2} ' description{myobs(1)} '/' description{myobs(2)}])
ylabel(['log_{2} ' description{myobs(3)} '/' description{myobs(4)}])

set(gca,'XLim',[-6 4],'YLim',[-2 4])

text(1,-1,['\chi^2/N = ' num2str(chi2/(2*length(erkaktratio(:))))])

