% Figure box: Explanation of fPCA
addpath('./Functions/')

close all

load('./Workspaces/harm_basis_50_to_end')
% load('./Workspaces/scores_puls')

% sites_all = [64];
sites_all = [4];
sigs = [1:60];
colmap = jet(length(sigs));
colind = sigs;
sorted_inds = sigs;
harmonicsRemoved = [1]; % specify harmonics that shall be removed from signal, e.g. for trend effects

% myscores = scores_puls(ismember(celltypes,sites_all),1);
myscores = edge_snr_score_pw_distdur(sites_all,[],0,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
myscores = myscores(sigs);
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
    load(['./Workspaces/site_' num2str(isite)])
    c_signal_single = log10(intensity);
    s = siteprop(isite);
    legstr{icount} = s.lig_name;
end


figure
setFigure(gcf,xfac,yfac,fontsize);

plot(harm_basis);legend(num2str((1:getnbasis(harm_basis))'));
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
[~, range_ind_min] = min(abs(timestamp - time_range(1) - 5));
[~, range_ind_max] = min(abs(timestamp - time_range(2) + 5));
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
for iplot = sigs(sorted_inds)
    legh = [legh plot(timestamp,c_signal_single(:,iplot),'o','Color',colmap(colind(iplot),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
    plot(times_fine_late,data_fpca_repr_fine(iplot,:),'Color',colmap(colind(iplot),:),'LineWidth',linewidth)
end


xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');
set(gca,'XLim',time_range)

figure
h = heatmap(data_fpca_repr_fine(sorted_inds,:), times_fine_late, [],-.015,.015,[],'Colormap','money','UseFigureColormap',false);
drawnow
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),50,'linear','extrap');
hold on;plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),200);
plot([xtick_new xtick_new],get(gca,'YLim'),':k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),500);
plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xlabel('time [min]')
ylabel('Localization');
set(gcf,'Position',[100 100 700 400]);
figure
setFigure(gcf,xfac,yfac,fontsize)

hold on

legh = [];
for isite = sigs(sorted_inds)
    legh = [legh plot(timestamp(range_ind),c_signal_woNharm(:,isite),'o-','Color',colmap(colind(isite),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
end


xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');
set(gca,'XLim',time_range)

figure
h2 = heatmap(c_signal_woNharm(:,sorted_inds)', timestamp(range_ind), [],-.015,.015,[],'Colormap','money','UseFigureColormap',false);

drawnow
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),50,'linear','extrap');
hold on;plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),200);
plot([xtick_new xtick_new],get(gca,'YLim'),':k');
xtick_new = interp1(str2num(char(get(gca,'XTickLabel'))),get(gca,'XTick'),500);
plot([xtick_new xtick_new],get(gca,'YLim'),'-k');
xlabel('time [min]')
ylabel('Pulsing');
set(gcf,'Position',[100 100 700 400]);
