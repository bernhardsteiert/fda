% Figure 1b (eyecatcher): Traces of BTC 100ng/ml in 184A1
addpath('./Functions/')

isite = 64;

load(['./Workspaces/site_' num2str(isite)])
c_signal = log10(intensity);

time_range = [50 510];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

figure;

s = siteprop(isite);

first_n = 7; % Plot first_n traces colored (2 times)

tmp_puls_strengths = edge_snr_score_pw_distdur(isite);
[tmp ind_tmp_pul_str] = sort(tmp_puls_strengths);
ind_isite = [ind_tmp_pul_str(end:-1:end-first_n+1) ind_tmp_pul_str(round(linspace(1,(length(ind_tmp_pul_str)-first_n),first_n)))];

plot(repmat(timestamp(range_ind),1,size(c_signal,2)),c_signal(range_ind,:),'g','color',[0.7 0.7 0.7])
hold on
plot(repmat(timestamp(range_ind),1,2*first_n),c_signal(range_ind,ind_isite))
% plot(timestamp(range_ind),nanmean(c_signal(range_ind,:),2),'color','k','LineWidth',2)
title([s.lig_name num2str(s.lig_dose) ' ng/ml'])

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

ylim = [-1 1]*.04;
plot([200 200],ylim,'k--')
text(125,.033,{'non-stationary','(deterministic)'},'HorizontalAlignment','center')
text(355,.033,{'stationary','(stochastic)'},'HorizontalAlignment','center')

set(gca,'XLim',time_range,'YLim',ylim)
set(gca,'XTick',50:50:500)