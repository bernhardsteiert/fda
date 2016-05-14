load('./Workspaces/harm_basis_130722_corrected_retracked_all_cleaned_late_newBTC')
load('./Workspaces/scores_early_5basis_noFGF_newBTC');
late_ws = load('./Workspaces/scores_puls_corrected_retracked_all_cleaned_newBTC');
close all;
myextension = '130722_corrected_retracked_all_cleaned';

sites_all = [17 57 64];
sigs = [26 35 99];
colind = [1 6 5];

colmap = [linspace(0,1,7)' ones(7,1) ones(7,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
linewidth = 0.75;

figure
hold on

c_signal_single = [];
scores_single = nan(5,length(sigs));

legstr = {};
for icount = 1:length(sites_all)
    isite = sites_all(icount);
    load(['./Workspaces/site_' num2str(isite) '_' myextension])
    c_signal_single(:,icount) = log10(intensity(:,sigs(icount)));
    scores_tmp = scores_early(:,celltypes==isite);
    scores_single(:,icount) = scores_tmp(:,sigs(icount));
    s = siteprop(isite);
    legstr{icount} = s.lig_name;
end

c_signal_single(isinf(c_signal_single)) = nan;
for i = 1:size(c_signal_single,2)
    if sum(isnan(c_signal_single(:,i))) > length(c_signal_single(:,i))-2
        c_signal_single(:,i) = 0;
    end
    c_signal_single(:,i) = interp1(timestamp(~isnan(c_signal_single(:,i))),c_signal_single(~isnan(c_signal_single(:,i)),i),timestamp);
    vec = ~isnan(c_signal_single(:,i))';
    rl = find(vec ~= [vec(2:end), vec(end)+1]);
    data =  vec(rl);
    rl(2:end) = rl(2:end) - rl(1:end-1);
    if ~data(1)
        c_signal_single(1:rl(1),i) = c_signal_single(rl(1)+1,i);
    end
    if ~data(end)
        c_signal_single(end-rl(end)+1:end,i) = c_signal_single(end-rl(end),i);
    end
end

legh = [];
for isite = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,isite),[markers{colind(isite)} '-'],'Color',colmap(colind(isite),:),'MarkerFaceColor','w','MarkerEdgeColor',colmap(colind(isite),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:')
end

legend(legh,legstr)

set(gca,'XLim',[50 700],'XTick',120:100:920,'XTickLabel',0:100:800)
% plot([200 200],get(gca,'YLim'),'k--')

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

time_range = [50.7 197.8];
times_fine = linspace(time_range(1),time_range(2),501);

times_fine(end+1) = 200;

axpos = get(gca,'Position');
xlim = get(gca,'XLim');
annotation('rectangle',[axpos(1:2) (times_fine(end)-xlim(1))./range(xlim) * axpos(3) axpos(4)])

% fPCA plot

sites_all_single = sites_all;
sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF

outliers = false(size(celltypes));
outlier_thres = .05;
nNeighbours = 8;
for iscore = 1:length(scores_early)
    sorted_dists = sort(sqrt(sum((repmat(scores_early(2:3,iscore),1,size(scores_early,2))-scores_early(2:3,:)).^2,1)));
    outliers(iscore) = sorted_dists(nNeighbours+1) > outlier_thres;
end


highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);

    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
end

resort = [2 3 4 1 6 5];

highdoses = highdoses(resort);

figure;

plot(scores_early(2,~ismember(celltypes,highdoses) & ~outliers),scores_early(3,~ismember(celltypes,highdoses) & ~outliers),'o','MarkerSize',1,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none')

hold on

color_ind = 1;
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
legstr = cell(length(sites_all_single),1);
legh = [];
markers = {'o','s','v','d','^','>'};
for isite = highdoses([6 2 3 4 5 1])
    sprop = siteprop(isite);
    
    scores = scores_early(:,celltypes == isite);
    if ismember(isite,sites_all_single)
        legstr{isite == sites_all_single} = sprop.lig_name(1:3);
        legh(isite == sites_all_single) = plot(scores_single(2,isite == sites_all_single),scores_single(3,isite == sites_all_single),markers{isite == highdoses},'MarkerSize',5,'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerEdgeColor','none');
    end
    plotEllipsis(scores(2,~outliers(celltypes == isite)),scores(3,~outliers(celltypes == isite)),colmap(isite == highdoses,:),.5);
end


xlim = [-.25 .35];
ylim = [-.05 .15];
set(gca,'XLim',xlim,'YLim',ylim)

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',1,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

ylabel('fPC3 score')
xlabel('fPC2 score')


% pulsing plot

puls_thres = 0.3;

% ylims = [[-.2 .1];[-.1 .15];[-.05 .15]];
figure
ylabel('Pulse Score')
hold on

for ilig = 1:length(sigs)
    tmp = repmat(scores_single(:,ilig),1,sum(celltypes==sites_all_single(ilig)));
    scores_puls = late_ws.scores_puls(late_ws.celltypes==sites_all_single(ilig));
    h = bar(ilig,scores_puls(sum(tmp == scores_early(:,celltypes==sites_all_single(ilig)),1)~= 0),1);
    set(h,'FaceColor',colmap(colind(ilig),:))
end

plot([.5 length(sigs)+.5],[1 1]*puls_thres,'k--')

set(gca,'XTick',1:size(late_ws.scores_puls,2),'XTickLabel',legstr)
% set(gca,'YLim',ylims(iplot,:))
