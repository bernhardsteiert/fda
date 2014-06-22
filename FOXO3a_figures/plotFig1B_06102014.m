close all;clear all;
load(['.\Workspaces\combined06102014']);
M = mydata(mydata.CellLine=='184A1' & mydata.LigandConcentration>1,:);
N{1} = mydata(mydata.CellLine=='184A1' & mydata.Time ==15 & mydata.LigandConcentration>1,:);
N{2} = mydata(mydata.CellLine=='184A1' & mydata.Time ==60 & mydata.LigandConcentration>1,:);
N{3} = mydata(mydata.CellLine=='184A1' & mydata.Time ==180 & mydata.LigandConcentration>1,:);
x = double(M.CoverNMean(M.Assay=='Parental',:));
y = double(M.CoverNMean(M.Assay=='Reporter',:));
RHO = corr(x,y,'type','Pearson');
figure;
mycolor = hsv(3);
for i=1:3
    x = double(N{i}.CoverNMean(N{i}.Assay=='Parental',:));
    y = double(N{i}.CoverNMean(N{i}.Assay=='Reporter',:));
    plot(x,y,'o','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor',mycolor(i,:)); hold on;
end
legend({'15MIN';'60MIN';'180MIN'})
%f = fit(x,y,'poly1');
%plot(f,x,y);
xlabel('Endogenous FoxO3a C/N Ratio');
ylabel('Reporter FoxO3a C/N Ratio');
title(['184A1 \rho =' num2str(RHO)]);


M = mydata(mydata.CellLine=='MCF10A' & mydata.LigandConcentration>1,:);
N{1} = mydata(mydata.CellLine=='MCF10A' & mydata.Time ==15 & mydata.LigandConcentration>1,:);
N{2} = mydata(mydata.CellLine=='MCF10A' & mydata.Time ==60 & mydata.LigandConcentration>1,:);
N{3} = mydata(mydata.CellLine=='MCF10A' & mydata.Time ==180 & mydata.LigandConcentration>1,:);
x = double(M.CoverNMean(M.Assay=='Parental',:));
y = double(M.CoverNMean(M.Assay=='Reporter',:));
RHO = corr(x,y,'type','Pearson');
figure;
mycolor = hsv(3);
for i=1:3
    x = double(N{i}.CoverNMean(N{i}.Assay=='Parental',:));
    y = double(N{i}.CoverNMean(N{i}.Assay=='Reporter',:));
    plot(x,y,'o','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor',mycolor(i,:)); hold on;
end
legend({'15MIN';'60MIN';'180MIN'})
%f = fit(x,y,'poly1');
%plot(f,x,y);
xlabel('Endogenous FoxO3a C/N Ratio');
ylabel('Reporter FoxO3a C/N Ratio');
title(['MCF10A \rho =' num2str(RHO)]);