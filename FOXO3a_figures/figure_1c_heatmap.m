% Figure box: Explanation of fPCA
addpath('./Functions/')

close all

load('./Workspaces/harm_basis_50_to_600')
% load('./Workspaces/scores_puls')
% myextension = '130722_corrected_retracked_all_cleaned';
myextension = '130722_corrected_retracked_all_paper_cleaned';

% sites_all = [64];
sites_all = 4;
nsigs = 40;

removeTraces = [4 7 13 18 23];
%removeTraces = [6 7 9 14 15 24 25 29 30 32 37 45 52];
sigs = setdiff(1:100,removeTraces);
sigs = sigs(1:nsigs);
colmap = jet(length(sigs));
colind = 1:nsigs;
harmonicsRemoved = [1]; % specify harmonics that shall be removed from signal, e.g. for trend effects

% myscores = scores_puls(ismember(celltypes,sites_all),1);
myscores = edge_snr_score_pw_distdur(sites_all,myextension,0,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
% myscores = myscores(sigs);
[tmp sorted_inds] = sort(myscores);
sorted_inds = sorted_inds(end:-1:1);

xfac = 1;
yfac = 1;
fontsize = 24;
linewidth = 2;

c_signal_single = [];
scores_single = nan(3,length(sigs));

for icount = 1:length(sites_all)
    isite = sites_all(icount);
    load(['./Workspaces/site_' num2str(isite) '_' myextension])
    c_signal_single = log10(intensity);
    s = siteprop(isite);
    legstr{icount} = s.lig_name;
end

nsigs = min([nsigs size(c_signal_single,2)]);
sigs = sigs(1:nsigs);
colind = 1:nsigs;

figure
setFigure(gcf,xfac,yfac,fontsize);

plot(harm_basis);legend(num2str((1:getnbasis(harm_basis))'));
set(gca,'XLim',[50 1000],'XTick',120:100:920,'XTickLabel',0:100:800)
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
time_range = getbasisrange(harm_basis);
range_ind_min = find((timestamp-time_range(1))>0);
range_ind_min = range_ind_min(1);
range_ind_max = find((timestamp-time_range(2))<0);
range_ind_max = range_ind_max(end);
range_ind = range_ind_min:range_ind_max;
times_fine_late = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201);
smoothed_additional = smooth_basis(timestamp(range_ind),c_signal_single(range_ind,:),harm_basis);
harm_eval = eval_basis(harm_basis,timestamp(range_ind));
harm_eval_fine = eval_basis(harm_basis,times_fine_late);
fitcoef = getcoef(smoothed_additional);

nharm = size(fitcoef,1);
remainingHarm = sort(setdiff(1:nharm,harmonicsRemoved));
data_fpca_repr_fine = fitcoef(remainingHarm,:)'*harm_eval_fine(:,remainingHarm)';
data_fpca_repr = fitcoef(harmonicsRemoved,:)'*harm_eval(:,harmonicsRemoved)';
c_signal_single(range_ind,:) = c_signal_single(range_ind,:)-data_fpca_repr';
data_fpca_repr = fitcoef(remainingHarm,:)'*harm_eval(:,remainingHarm)';
c_signal_woNharm = c_signal_single(range_ind,:)-data_fpca_repr';

figure
setFigure(gcf,xfac,yfac,fontsize)

hold on
legh = [];
for iplot = 1:nsigs
    legh = [legh plot(timestamp,c_signal_single(:,sigs(iplot)),'o','Color',colmap(colind(iplot),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
    plot(times_fine_late,data_fpca_repr_fine(sigs(iplot),:),'Color',colmap(colind(iplot),:),'LineWidth',linewidth)
end


xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');
set(gca,'XLim',[50 600],'XTick',120:100:520,'XTickLabel',0:100:400)

figure
colmap2 = [];
for i = 1:250
    colmap2 = [colmap2; hsv2rgb([0 1-(i-1)/250 1])];
end
colmap2 = [colmap2; [1 1 1]];
for i = 1:250
    colmap2 = [colmap2; hsv2rgb([2/3 i/250 1])];
end
h = heatmap(data_fpca_repr_fine(sorted_inds(sigs),:),[],[],[],'Colormap',colmap2,'UseFigureColormap',false);
set(get(h,'Parent'),'CLim',[-.015 .015])
xtickfac = size(data_fpca_repr_fine,2)/(timestamp(range_ind(end))-timestamp(range_ind(1)));
xticklab = 70:50:600;
set(gca,'XLim',[0 size(data_fpca_repr_fine,2)]+.5,'XTick',(xticklab-timestamp(range_ind(1)))*xtickfac,'XTickLabel',xticklab-120)
hold on
% xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),-70,'linear','extrap');
% plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),0);
plot([xtick_new xtick_new],get(gca,'YLim'),':k');
% xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),480);
% plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xlabel('time [min]')
ylabel('Localization');
set(gcf,'Position',[100 100 400 250]);

% Old code to adjust heatmap x-ticks ... not working atm
% xtick = round(linspace(1,length(timestamp),11));
% xticklab_new = [50 120 200:100:700];
% xticklab = timestamp(xtick);
% xtick_new = interp1(timestamp(xtick),xtick,xticklab_new);
% set(get(h2,'Parent'),'XTick',xtick_new)
% set(get(h2,'Parent'),'XTickLabel',num2str(xticklab_new'))

figure
setFigure(gcf,xfac,yfac,fontsize)

hold on

legh = [];
for isite = 1:nsigs
    legh = [legh plot(timestamp(range_ind),c_signal_woNharm(:,sigs(isite)),'o-','Color',colmap(colind(isite),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
end


xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');
set(gca,'XLim',[50 600],'XTick',120:100:520,'XTickLabel',0:100:400)

figure
% h2 = heatmap(c_signal_woNharm(:,sorted_inds(sigs))', timestamp(range_ind), [],-.015,.015,[],'Colormap',colmap2,'UseFigureColormap',false);
h2 = heatmap(c_signal_woNharm(:,sorted_inds(sigs))',[],[],[],'Colormap',colmap2,'UseFigureColormap',false);
set(get(h2,'Parent'),'CLim',[-.015 .015])
xtickfac = size(c_signal_woNharm,1)/(timestamp(range_ind(end))-timestamp(range_ind(1)));
xticklab = 70:50:600;
set(gca,'XLim',[0 size(c_signal_woNharm,1)]+.5,'XTick',(xticklab-timestamp(range_ind(1)))*xtickfac,'XTickLabel',xticklab-120)

drawnow
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),-70,'linear','extrap');
hold on;plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),0);
plot([xtick_new xtick_new],get(gca,'YLim'),':k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),480);
plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xlabel('time [min]')
ylabel('Pulsing');
set(gcf,'Position',[100 100 400 250]);
