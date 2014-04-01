% Figure box: Explanation of fPCA
addpath('./Functions/')

close all

load('./Workspaces/harm_basis_50_to_end')
% load('./Workspaces/scores_puls')

sites_all = [64];
sigs = [1:50];
colmap = jet(length(sigs));
colind = sigs;
sorted_inds = sigs;

% myscores = scores_puls(ismember(celltypes,sites_all),1);
% myscores = myscores(sigs);
% [tmp sorted_inds] = sort(myscores);
% sorted_inds = sorted_inds(end:-1:1);

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
setFigure(gcf,xfac,yfac,fontsize)
plot(harm_basis)


time_range = getbasisrange(harm_basis);
[tmp range_ind_min] = min(abs(timestamp - time_range(1) - 5));
[tmp range_ind_max] = min(abs(timestamp - time_range(2) + 5));
range_ind = range_ind_min:range_ind_max;
times_fine_late = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201);
smoothed_additional = smooth_basis(timestamp(range_ind),c_signal_single(range_ind,:),harm_basis);
harm_eval = eval_basis(harm_basis,timestamp(range_ind));
harm_eval_fine = eval_basis(harm_basis,times_fine_late);
fitcoef = getcoef(smoothed_additional);
data_fpca_repr = fitcoef'*harm_eval';
data_fpca_repr_fine = fitcoef'*harm_eval_fine';
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

figure
h = heatmap(data_fpca_repr_fine(sorted_inds,:), times_fine_late, [],-.03,.025,[],'Colormap','money','UseFigureColormap',false);
drawnow
xtick = round(linspace(1,length(times_fine_late),11));
xticklab_new = [50 120 200:100:700];
xticklab = times_fine_late(xtick);
xtick_new = interp1(times_fine_late(xtick),xtick,xticklab_new);
set(get(h,'Parent'),'XTick',xtick_new)
set(get(h,'Parent'),'XTickLabel',num2str(xticklab_new'))

xlabel('time [min]')
ylabel('Localization');

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

figure
h2 = heatmap(c_signal_woNharm(:,sorted_inds)', timestamp, [],-.01,.01,[],'Colormap','money','UseFigureColormap',false);
drawnow
xtick = round(linspace(1,length(timestamp),11));
xticklab_new = [50 120 200:100:700];
xticklab = timestamp(xtick);
xtick_new = interp1(timestamp(xtick),xtick,xticklab_new);
set(get(h2,'Parent'),'XTick',xtick_new)
set(get(h2,'Parent'),'XTickLabel',num2str(xticklab_new'))

xlabel('time [min]')
ylabel('Pulsing');

