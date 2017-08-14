% Figure 2a: Linear correlation of pERK/pAKT vs. pS294/pS253
addpath('./Functions/')

% [westerndata,description,raw] = xlsread('./Workspaces/westernBlotData_final2');
[westerndata,description,raw] = xlsread('./Workspaces/westernBlotData_final3');
experiments_investigated = 1:2;

liglabels = {'NS','EGF','IGF','HRG','HGF','EPR','BTC'};
tmpligs = [1 2 3 5 6 7 8];
exind = 1;
timeind = 2;
ligind = 3;
obsind = 4:7;
myobs = [4 5 6 7];
normind = 8;
if all(isnan(westerndata(:,ligind)))
    for i = 1:length(liglabels)
        ind = strcmp(liglabels{i},raw(2:end,ligind));
        westerndata(ind,ligind) = tmpligs(i);
    end
end
uni_ligs = unique(westerndata(:,ligind));
resort = [1 3 4 5 2 7 6];
uni_ligs = uni_ligs(resort);
liglabels = liglabels(resort);

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

% Scaling of experiments
scaledLigs = data(data(:,exind) == experiments_investigated(1),ligind);
scaledTime = data(data(:,exind) == experiments_investigated(1),timeind);
pinit = data(data(:,exind) == experiments_investigated(1),myobs);
pinit = [ones(length(experiments_investigated)-1,length(myobs)); pinit];
[optRes,~,~,~,~,~,jacobian] = lsqnonlin(@(p) res(p,data(:,[exind timeind myobs])),pinit);

Cov=inv(jacobian'*jacobian);
scaledData = optRes(length(experiments_investigated):end,:);
scaledStd = sqrt(diag(Cov));
scaledStd = reshape(scaledStd,size(optRes));
scaledStd = scaledStd(length(experiments_investigated):end,:);

rowstocols = 0.5;
nrows = ceil(size(scaledData,2)^rowstocols);
ncols = ceil(size(scaledData,2) / nrows);

% Plotting time-courses
colmap = [linspace(0,1,length(uni_ligs))' ones(length(uni_ligs),1) ones(length(uni_ligs),1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
colmap = [colmap; [0 0 0]];

figure
legh = [];
stdwidth = 1;
for iplot = 1:size(scaledData,2)
    subplot(nrows,ncols,iplot)
    hold on
    colcount = 1;
    for ilig = uni_ligs(2:end)'
        myind = logical([1; scaledLigs(2:end) == ilig]); % Plot common timepoint zero for all ligands
        if sum(~isnan(scaledData(myind,iplot))) > 1 % not only NaN
            tmpx = [scaledTime(myind); flipud(scaledTime(myind))];
            tmpy = [scaledData(myind,iplot) + scaledStd(myind,iplot)./stdwidth; flipud(scaledData(myind,iplot) - scaledStd(myind,iplot)./stdwidth)];
            
            ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
            set(ltmp, 'FaceColor', colmap(colcount,:)*0.1+0.9, 'EdgeColor', colmap(colcount,:)*0.1+0.9);
            ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
            set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(colcount,:)*0.3+0.7);
                                    
                                    
            legh = [legh plot(scaledTime(myind),scaledData(myind,iplot),'x-','Color',colmap(colcount,:))];
            colcount = colcount + 1;
        end
    end
    title(description{myobs(iplot)})
    ylabel('log_2 fold change [au]')
    xlabel('time [min]')
    set(gca,'XLim',[-10 490])
     set(gca,'YLim',[-7 1],'YTick',[-6:2:0])
end

legend(legh,liglabels{2:end})

% Plotting correlation diagram
markers = {'o','s','v','d','^','>','<'};
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
        legh = [legh errorbar(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4),sqrt(sum(scaledStd(myind,3:4).^2./stdwidth,2)),markers{colcount},'MarkerFaceColor',colmap(colcount,:),'MarkerEdgeColor',colmap(colcount,:),'Color',colmap(colcount,:))];
        h = herrorbar(scaledData(myind,1)-scaledData(myind,2),scaledData(myind,3)-scaledData(myind,4),sqrt(sum(scaledStd(myind,1:2).^2./stdwidth,2)));
        set(h,'Color',colmap(colcount,:))
        set(h(2),'LineStyle', 'none')
        colcount = colcount + 1;
    end
end
axb = lsqnonlin(@(axb) [(erkaktratio(:)-(foxositesratio(:)-axb(2))./axb(1))./(erkaktratio_std(:)); (foxositesratio(:)-axb(1)*foxositesratio(:)-axb(2))./(foxositesratio_std(:))],[1 0]);
plot([min(min(erkaktratio)) max(max(erkaktratio))],[min(min(erkaktratio)) max(max(erkaktratio))]*axb(1) + axb(2),'k--','LineWidth',2)
chi2 = sum(([(erkaktratio(:)-(foxositesratio(:)-axb(2))./axb(1))./(erkaktratio_std(:)); (foxositesratio(:)-axb(1)*foxositesratio(:)-axb(2))./(foxositesratio_std(:))]).^2);

legend(legh,liglabels{2:end},'Location','NorthWest')
xlabel(['log_{2} ' description{1,myobs(1)} '/' description{1,myobs(2)}])
ylabel(['log_{2} ' description{1,myobs(3)} '/' description{1,myobs(4)}])

set(gca,'XLim',[-7.5 5],'YLim',[-3 3.5],'XTick',[-6:2:4])

text(1,-1,['\chi^2/N = ' num2str(chi2/(2*length(erkaktratio(:))))])

return
% Write xls file
description_scaled = {'Time','LigInd','pERK','pAKT','pS294','pS253','pERK_std','pAKT_std','pS294_std','pS253_std'};
scaled_final = [scaledTime scaledLigs scaledData full(scaledStd)];
fprintf('%s,',description_scaled{:})
fprintf('\n')
dlmwrite('Workspaces/westernBlot_scaled.csv',scaled_final)