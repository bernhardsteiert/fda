% Figure 2b: Traces of all ligands at 100ng/ml in 184A1
addpath('./Functions/')

% sites_all = [17 37 44 4 57 64];
sites_all = [17 37 44 64 57 4];
sites_akti = [19 39 42 62 59 2];

times = cell(0);
signals = cell(0);
celltype = [];
legstr = {};

for isite = [sites_all sites_akti]
    load(['./Workspaces/site_' num2str(isite)])
    times{end+1} = timestamp;  
    signals{end+1} = log10(intensity);
    celltype = [celltype ones(1,size(intensity,2))*isite];
    s = siteprop(isite);
    legstr{end+1} = s.lig_name;
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);

time_range = [50 605];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max+1;

figure

ncols = 2;
nrows = 1;

colmap = [linspace(0,1,length(sites_all)+1)' ones(length(sites_all)+1,1) ones(length(sites_all)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));

ylim = [-1 1]*.05;

subplot(nrows,ncols,1)
hold on

resort = [1 2 3 6 4 5];
legh = [];
for iplot = 1:length(sites_all)
    isite = sites_all(resort(iplot));
    legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == isite),2),'Color',colmap(iplot,:))];
end
title('No Drug')
set(gca,'XLim',time_range,'YLim',ylim/1.5)
set(gca,'XTick',120:100:520,'XTickLabel',0:100:400,'YTick',-.03:.01:.03)
xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]')

subplot(nrows,ncols,2)
hold on

resort = [1 2 3 6 4 5];
legh = [];
for iplot = 1:length(sites_akti)
    isite = sites_akti(resort(iplot));
    legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == isite),2),'Color',colmap(iplot,:))];
end
title('AKTi')
set(gca,'XLim',time_range,'YLim',ylim/1.5)
set(gca,'XTick',120:100:520,'XTickLabel',0:100:400,'YTick',-.03:.01:.03)
legend(legh,legstr{resort})